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
using GridapGmsh
using LinearAlgebra

to = TimerOutput()
@timeit to "variables" begin
    # Variables
    g = 9.81        # [kg/s²] gravitational constant
    ϕ₀(x) = 1.0 + im * 0.0
end

function setup_domain(B,D,pmid)
    geo = quadrilateral(;x0=pmid,d1=VectorValue(B/2,0.0),d2=VectorValue(0.0,D))
    # model = CartesianDiscreteModel(pmin, pmax, partition)
    model = GmshDiscreteModel("data/meshes/background6.msh")
    # backgroundmesh5.msh for high res first order elements
    model,geo
end # function

function cutting_model(model,geo)
    @timeit to "embedding" begin
        cutgeo = cut(model, !geo)
        cutgeo_facets = cut_facets(model, !geo)
    end
    cutgeo, cutgeo_facets
end # function

# List of H/R ratios corresponding to each mesh model in the `models` array
#HR_ratios = ["0", "0.342", "0.643", "0.809", "0.906", "0", "-0.259", "-0.643", "-0.809", "-0.906"]
HR_ratios = [0.0]
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
    n₁ = 100
    n₂ = 25
    D = 0.1
    B = 2*D
    h = 40*D       # [m] depth
    # ω = √(0.1*2*g/B)
    k_max = 1.0
    λ = 2*π/k_max 
    L₁ = 2*λ
    for HR in HR_ratios
        model, geo = setup_domain(B,D,Point(0.0,-D))
        push!(models, model)
        push!(geos, geo)
    end # for
end

function setup_spaces(order::Int64, model::DiscreteModel, Ω, cutgeo)
    reffeᵩ = ReferenceFE(lagrangian, Float64, order)
    Wstd = FESpace(Ω, reffeᵩ, vector_type=Vector{ComplexF64})
    threshold = 1.0
    strategy = AggregateCutCellsByThreshold(threshold)
    aggregates = aggregate(strategy, cutgeo, cutgeo.geo)
    W = AgFEMSpace(Wstd,aggregates)
    Φ = TrialFESpace(W)
    V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{2, ComplexF64})
    U = TrialFESpace(V)
    W,Φ,V,U
end # function

function setup_spaces(order::Int64, model::DiscreteModel, Ω)
    reffeᵩ = ReferenceFE(lagrangian, Float64, order)
    W = FESpace(Ω, reffeᵩ, vector_type=Vector{ComplexF64})
    Φ = TrialFESpace(W)
    V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{2, ComplexF64})
    U = TrialFESpace(V)
    W,Φ,V,U
end # function

function setup_interiors(model,cutgeo,cutgeo_facets,degree)
    Ω = Interior(cutgeo, PHYSICAL)
    Ω⁻act = Interior(cutgeo, ACTIVE)
    Γ = EmbeddedBoundary(cutgeo)
    nΓ = get_normal_vector(Γ)
    Γf⁻ = BoundaryTriangulation(cutgeo_facets, tags=["surface"])
    dΓf⁻ = Measure(Γf⁻, degree)
    dΩ = Measure(Ω, degree)
    dΓ = Measure(Γ, degree)
    Γi = BoundaryTriangulation(model, tags=["inlet"])
    dΓi = Measure(Γi, degree)
    Γo = BoundaryTriangulation(model, tags=["outlet"])
    dΓo = Measure(Γo, degree)
    E⁰ = GhostSkeleton(cutgeo)
    dE⁰ = Measure(E⁰,degree)
    nE⁰ = get_normal_vector(E⁰)
    (Ω,Ω⁻act),(Γ,Γf⁻,Γi,Γo,E⁰),(nΓ,nE⁰),(dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰)
end # function

function setup_interiors(model,cutgeo,degree;method="SBM")
    Ω = Interior(cutgeo, IN)
    writevtk(Ω,"model")
    Γ = Interface(Interior(cutgeo,ACTIVE_OUT),Ω).⁻
    dΩ = Measure(Ω, degree)
    dΓ = Measure(Γ, degree)   
    nΓ = get_normal_vector(Γ)
    Γf⁻ = BoundaryTriangulation(Ω, tags=["surface"])
    dΓf⁻ = Measure(Γf⁻, degree)
    Γi = BoundaryTriangulation(model, tags=["inlet"])
    dΓi = Measure(Γi, degree)
    Γo = BoundaryTriangulation(model, tags=["outlet"])
    dΓo = Measure(Γo, degree)
    (Ω),(Γ,Γf⁻,Γi,Γo),(nΓ),(dΩ,dΓ,dΓf⁻,dΓi,dΓo)
end # function

function d(x,t)
    dx = maximum([-0.1-x[1], 0.0, x[1]-0.1])
    dy = maximum([-0.1-x[2], 0.0, x[2]-0.1])
    VectorValue(-dx,dy)
end # function

function n(x,t)
    dx = maximum([-0.1-x[1], 0.0, x[1]-0.1])
    dy = maximum([-0.1-x[2], 0.0, x[2]-0.1])
    dist = √(dx^2+dy^2)
    VectorValue(-dx/dist,dy/dist)
end # function

d(t) = x -> d(x,t)
n(t) = x -> n(x,t)

function main_agfem(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,nΓ,W,Φ,V,U)


    # Matrix assemblies #agfem + cutfem
    a_wϕ(ϕ, w) = ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*(w*ϕ) )dΓf⁻ - ∫( im*k*ϕ*w )dΓo 
    a_vϕ(ϕ, v) = ∫( (v⋅nΓ) * (  im*ω*ϕ )* (-1.0))dΓ
    a_wu(u, w) =  ∫(w * im*ω*(u⋅nΓ) )dΓ

    A_wϕ = assemble_matrix(a_wϕ, Φ, W)
    A_wu = assemble_matrix(a_wu, U, W)
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)
    # A_vu = assemble_matrix(a_vu, U, V)

    AA = A_vϕ * (inv(Matrix(A_wϕ))) * A_wu
    ρV = (12*4-∑(∫(1.0)dΩ))#D*B/2
    Mₐ = (real(AA)  )/ ρV / (ω^2) 
    Cₐ = imag(AA) / (ω^2) / ρV
    return Mₐ, Cₐ
end # function

function main_cutfem(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰,nΓ,nE⁰,W,Φ,V,U)
    # Matrix assemblies #agfem + cutfem
    a_wϕ(ϕ, w) = ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*(w*ϕ) )dΓf⁻ - ∫( im*k*ϕ*w )dΓo #+∫( 0.1 * jump(∇(w)⋅nE⁰) * jump( ∇(ϕ)⋅nE⁰) )dE⁰
    a_vϕ(ϕ, v) = ∫( (v⋅nΓ) * (  im*ω*ϕ )* (-1.0))dΓ
    a_wu(u, w) =  ∫(w * im*ω*(u⋅nΓ) )dΓ

    A_wϕ = assemble_matrix(a_wϕ, Φ, W)
    A_wu = assemble_matrix(a_wu, U, W)
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)
    # A_vu = assemble_matrix(a_vu, U, V)

    AA = A_vϕ * (inv(Matrix(A_wϕ))) * A_wu
    ρV = (12*4-∑(∫(1.0)dΩ))#D*B/2
    Mₐ = (real(AA)  )/ ρV / (ω^2) 
    Cₐ = imag(AA) / (ω^2) / ρV
    return Mₐ, Cₐ
end # function

function main_sbm(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,nΓ,W,Φ,V,U )
    a_wϕ(ϕ, w) = ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*(w*ϕ) )dΓf⁻ - ∫( im*k*ϕ*w )dΓo + ∫(w*(n(0)⋅nΓ)*((∇∇(ϕ)⋅d(0) + ∇(ϕ))⋅n(0)) - w* ∇(ϕ)⋅nΓ)dΓ # +∫( 0.1 * jump(∇(w)⋅nE⁰) * jump( ∇(ϕ)⋅nE⁰) )dE⁰
    a_vϕ(ϕ, v) = ∫( (v⋅n(0)) *(n(0)⋅nΓ) * (  im*ω*ϕ + im*ω*∇(ϕ)⋅d(0) )* (-1.0))dΓ
    a_wu(u, w) =  ∫(w * im*ω*(u⋅n(0))*(n(0)⋅nΓ) )dΓ

    A_wϕ = assemble_matrix(a_wϕ, Φ, W)
    A_wu = assemble_matrix(a_wu, U, W)
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)
    # A_vu = assemble_matrix(a_vu, U, V)

    # @time AA = A_vϕ / cholA_wϕ * A_wu
    AA = A_vϕ * (inv(Matrix(A_wϕ))) * A_wu
    ρV = (12*4-∑(∫(1.0)dΩ))#/0.01 *D*B/2
    # @show (12*4-∑(∫(1.0)dΩ))
    # @show (12*4-∑(∫(1.0)dΩ))/0.01
    # @show ∑(∫(1.0)dΩ)
    # @show ∑(∫(1.0)dΓ)/0.2
    # @show 0.2/∑(∫(1.0)dΓ)
    # @show ∑(∫(1.0)dΓ)
    # stop
    # ρV = D*B/2
    Mₐ = (real(AA)  )/ ρV / (ω^2) 
    Cₐ = imag(AA) / (ω^2) / ρV
    return Mₐ, Cₐ
end # function

KRs = 0.1:0.1:2.5
mass_plots = plot(title="Added Mass vs KR", xlabel="KR [-]", ylabel="Added Mass Mₐ [-]")
damping_plots = plot(title="Added Damping vs KR", xlabel="KR [-]", ylabel="Added Damping Cₐ [-]")

# Loop over each model and compute added mass and damping
for (i, model) in enumerate(models)
    added_mass = []
    added_damping = []
    geo = geos[i]
    cutgeo, cutgeo_facets = cutting_model(model,geo)
    order = 2
    degree = 2 * order

    # agfem
    # (Ω,Ω⁻act),(Γ,Γf⁻,Γi,Γo,E⁰),(nΓ,nE⁰),(dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰) = setup_interiors(model,cutgeo,cutgeo_facets,degree)
    # W,Φ,V,U = setup_spaces(order, model, Ω⁻act, cutgeo)

    # cutfem
    # (Ω,Ω⁻act),(Γ,Γf⁻,Γi,Γo,E⁰),(nΓ,nE⁰),(dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰) = setup_interiors(model,cutgeo,cutgeo_facets,degree)
    # W,Φ,V,U = setup_spaces(order, model, Ω⁻act)

    # sbm
    (Ω),(Γ,Γf⁻,Γi,Γo),(nΓ),(dΩ,dΓ,dΓf⁻,dΓi,dΓo) = setup_interiors(model,cutgeo,degree)
    W,Φ,V,U = setup_spaces(order, model, Ω)

    for KR in KRs
        k = KR/D
        ω = √(k * g)
        # Mₐ, Cₐ = main_agfem(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,nΓ,W,Φ,V,U)
        # Mₐ, Cₐ = main_cutfem(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰,nΓ,nE⁰,W,Φ,V,U)
        Mₐ, Cₐ = main_sbm(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,nΓ,W,Φ,V,U)
        push!(added_mass, Mₐ[2, 2])  # Extract relevant component for plotting
        push!(added_damping, Cₐ[2, 2])
        println(Mₐ[2, 2],",",Cₐ[2, 2])
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