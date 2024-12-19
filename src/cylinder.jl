module CONFORMAL_2D_CYLINDER_SPECTRAL
using Gridap
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.Fields
using GridapEmbedded.Interfaces
using WriteVTK
using GridapEmbedded
using GridapEmbedded.LevelSetCutters
using Gridap.Arrays
using Gridap.Adaptivity
using Gridap.ODEs
using TimerOutputs
using DrWatson
using Plots
using Gridap.FESpaces
# using GridapGmsh

to = TimerOutput()
@timeit to "variables" begin
    # Variables
    g = 9.81        # [kg/s²] gravitational constant
    ϕ₀(x) = 1.0 + im * 0.0
end

function setup_domain(L₁,h,partition,R,pmid)
    pmin = Point(-L₁/2, -h)
    pmax = Point(L₁/2, 0.0)
    geo = disk(R, x0=pmid)
    model = CartesianDiscreteModel(pmin, pmax, partition)
    labels_Γ = get_face_labeling(model)  # get the face labeling of model_Γ 
    add_tag_from_tags!(labels_Γ, "seabed", [1,2,5])  
    add_tag_from_tags!(labels_Γ, "inlet", [3, 7])
    add_tag_from_tags!(labels_Γ, "outlet", [4, 8]) 
    add_tag_from_tags!(labels_Γ, "surface", [6]) 
    model,geo
end # function

# List of H/R ratios corresponding to each mesh model in the `models` array
#HR_ratios = ["0", "0.342", "0.643", "0.809", "0.906", "0", "-0.259", "-0.643", "-0.809", "-0.906"]
HR_ratios = [0, 0.342, 0.643, 0.809, 0.906]
# HR_ratios = [0, 0.342, 0.643, 0.809, 0.906, -0.259, -0.643, -0.809, -0.906]

@timeit to "domain" begin
    # # Load mesh models for each H/R ratio
    # models = [
    #     GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_0_R_0_1.msh"),
    #     GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_0_342_R_0_1.msh"),
    #     GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_0_643_R_0_1.msh"),
    #     GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_0_809_R_0_1.msh"),
    #     GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_0_906_R_0_1.msh"),
    #     #GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_0_R_0_1.msh"),
    #     #GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_negative_0_259.msh"),
    #     #GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_negative_0_643.msh"),
    #     #GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_negative_0_809.msh"),
    #     #GmshDiscreteModel("Cylinder/data/meshes/HR_ratio/cylinder_HR_negative_0_906.msh")
    # ]
    models = []
    geos = []
    n₁ = 80
    n₂ = 40
    L₁ = 0.8
    h = 0.4
    R = 0.1
    for HR in HR_ratios
        model, geo = setup_domain(L₁,h,(n₁,n₂),R,Point(0.0,-R*HR))
        push!(models, model)
        push!(geos, geo)
    end # for
end

function main(KR, model, geo)
    @timeit to "embedding" begin
        cutgeo = cut(model, !geo)
        cutgeo_facets = cut_facets(model, !geo)
        # Ω = Interior(model) 
        Ω = Interior(cutgeo, PHYSICAL)
        Ω⁻act = Interior(cutgeo, ACTIVE)
        Γ = EmbeddedBoundary(cutgeo)
        nΓ = get_normal_vector(Γ)

        order = 1
        degree = 2 * order
        Γf⁻ = BoundaryTriangulation(cutgeo_facets, tags=["surface"])
        dΓf⁻ = Measure(Γf⁻, degree)
        nΓf = get_normal_vector(Γf⁻)
        dΩ = Measure(Ω, degree)
        dΓ = Measure(Γ, degree)

        Γsb = BoundaryTriangulation(model, tags=["seabed"])
        dΓsb = Measure(Γsb, degree)
        nΓsb = get_normal_vector(Γsb)
        Γi = BoundaryTriangulation(model, tags=["inlet"])
        dΓi = Measure(Γi, degree)
        nΓi = get_normal_vector(Γi)
        Γo = BoundaryTriangulation(model, tags=["outlet"])
        dΓo = Measure(Γo, degree)
        nΓo = get_normal_vector(Γo)
    end

    @timeit to "spaces" begin
        # Define FE spaces
        reffeᵩ = ReferenceFE(lagrangian, Float64, order)
        # reffeu = ReferenceFE(lagrangian, VectorValue{2, Float64}, order)
        Wstd = FESpace(Ω⁻act, reffeᵩ, vector_type=Vector{ComplexF64})
        Φstd = TrialFESpace(Wstd)
        threshold = 1.0
        strategy = AggregateCutCellsByThreshold(threshold)
        aggregates = aggregate(strategy, cutgeo, geo, OUT)
        W = AgFEMSpace(Wstd,aggregates)
        Φ = TrialFESpace(W)
        # E = FESpace(Ω⁻act, reffeᵩ, vector_type=Vector{ComplexF64})
        # Η = TrialFESpace(E)
        V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{2, ComplexF64})
        U = TrialFESpace(V)

        # X = MultiFieldFESpace([Φ, U])
        # Y = MultiFieldFESpace([W, V])
    end

    k = KR / R
    ω = √(k * g)
    # nz = VectorValue(0.0, 1.0)
    # m = 0.96 + im * 0.0
    # η₀ = 0.05
    # ∇ₙϕd(x) = -im * ω * η₀ * exp(im * k * x[2])

    # Matrix assemblies
    a_wϕ(ϕ, w) = ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*(w*ϕ) )dΓf⁻ - ∫( im*k*ϕ*w )dΓi - ∫( im*k*ϕ*w )dΓo
    a_vϕ(ϕ, v) = ∫( (v⋅nΓ) * (  im*ω*ϕ )* (-1.0))dΓ
    a_wu(u, w) =  ∫(w * im*ω*(u⋅nΓ) )dΓ
    # a_vu(u, v) = ∫(v ⋅ ( g*(u⋅nz)*nz + m*(ω^2)*u))dΓ

    A_wϕ = assemble_matrix(a_wϕ, Φ, W)
    A_wu = assemble_matrix(a_wu, U, W)
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)
    # A_vu = assemble_matrix(a_vu, U, V)

    AA = A_vϕ * (inv(Matrix(A_wϕ))) * A_wu
    ρV = (π * R^2)
    Mₐ = (real(AA)  )/ ρV / (ω^2) 
    Cₐ = imag(AA) / (ω^2) / ρV
    return Mₐ, Cₐ
end

KRs = 0.05:0.025:2.0
mass_plots = plot(title="Added Mass vs KR", xlabel="KR", ylabel="Added Mass Mₐ")
damping_plots = plot(title="Added Damping vs KR", xlabel="KR", ylabel="Added Damping Cₐ")

# Loop over each model and compute added mass and damping
for (i, model) in enumerate(models)
    added_mass = []
    added_damping = []
    geo = geos[i]
    for KR in KRs
        Mₐ, Cₐ = main(KR, model, geo)
        push!(added_mass, Mₐ[2, 2])  # Extract relevant component for plotting
        push!(added_damping, Cₐ[2, 2])
    end

    # Plot added mass and damping for the current model with H/R ratio label
    plot!(mass_plots, KRs, added_mass, label="H/R = $(HR_ratios[i])")
    plot!(damping_plots, KRs, added_damping, label="H/R = $(HR_ratios[i])")
end

# Display and save the plots
display(mass_plots)
display(damping_plots)
# savefig(mass_plots, "Added_Mass_vs_KR.png")
# savefig(damping_plots, "Added_Damping_vs_KR.png")

end