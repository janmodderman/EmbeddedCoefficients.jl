module TRIANGLE
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
    R = 0.1
    ϕ₀(x) = 1.0 + im * 0.0
end

function setup_domain(D,pmid)
    φ=30*π/180
    geo = quadrilateral(;x0=pmid,d1=VectorValue(D*tan(φ),D),d2=VectorValue(0.0,D))
    model = GmshDiscreteModel("data/meshes/background_shapes.msh")
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
    models = []
    geos = []
    for HR in HR_ratios
        model, geo = setup_domain(R,Point(0.0,-R))
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
    writevtk(Ω,"modelagfem")
    writevtk(Interior(model),"background")
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
    writevtk(Ω,"modelsbm")
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

function d(x,t,pmid)
    # dist=√((x[1]-pmid[1])^2+(x[2]-pmid[2])^2)-R
    dr = pmid-x
    dist = abs(√(dr[1]^2+dr[2]^2)-R)
    dist.*n(x,t,pmid)
end # function

function n(x,t,pmid)
    # dist=√((x[1]-pmid[1])^2+(x[2]-pmid[2])^2)
    dr = pmid-x
    dist = √(dr[1]^2+dr[2]^2)
    dr./dist
end # function

d(t,pmid) = x -> d(x,t,pmid)
n(t,pmid) = x -> n(x,t,pmid)

# @show n(Point(0.0,-0.1), 0, Point(0.0,-0.906*0.1))
# @show n(Point(0.1,-0.906*0.1), 0, Point(0.0,-0.906*0.1))
# @show n(Point(0.1,-0.1), 0, Point(0.0,-0.906*0.1))

# @show d(Point(0.0,-0.1), 0, Point(0.0,-0.906*0.1))
# @show d(Point(0.1,-0.906*0.1), 0, Point(0.0,-0.906*0.1))
# @show d(Point(0.1,-0.1), 0, Point(0.0,-0.906*0.1))
# stop

function main_agfem(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,nΓ,W,Φ,V,U)


    # Matrix assemblies #agfem + cutfem
    a_wϕ(ϕ, w) = ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*(w*ϕ) )dΓf⁻ - ∫( im*k*ϕ*w )dΓo 
    a_vϕ(ϕ, v) = ∫( (v⋅nΓ) * (  im*ω*ϕ )* (-1.0))dΓ
    a_wu(u, w) =  ∫(w * im*ω*(u⋅nΓ) )dΓ
    # rhs(w) = ∫(w*1.0)dΩ 

    A_wϕ = assemble_matrix(a_wϕ, Φ, W)
    A_wu = assemble_matrix(a_wu, U, W)
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)
    # A_vu = assemble_matrix(a_vu, U, V)
    # onesrhs=assemble_vector(rhs,W)

    AA = A_vϕ * (inv(Matrix(A_wϕ))) * A_wu
    # AA = A_vϕ * (A_wϕ\onesrhs) * A_wu
    ρV = R^2*tan(30*π/180)/2
    Mₐ = (real(AA)  )/ ρV / (ω^2) 
    Cₐ = imag(AA) / (ω^2) / ρV
    return Mₐ, Cₐ
end # function

function main_cutfem(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰,nΓ,nE⁰,W,Φ,V,U)
    # Matrix assemblies #agfem + cutfem
    a_wϕ(ϕ, w) = ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*(w*ϕ) )dΓf⁻ - ∫( im*k*ϕ*w )dΓo + ∫( 0.1*(0.0007)^(2*1+1) * jump(∇(w)⋅nE⁰) * jump( ∇(ϕ)⋅nE⁰) )dE⁰
    a_vϕ(ϕ, v) = ∫( (v⋅nΓ) * (  im*ω*ϕ )* (-1.0))dΓ
    a_wu(u, w) =  ∫(w * im*ω*(u⋅nΓ) )dΓ

    A_wϕ = assemble_matrix(a_wϕ, Φ, W)
    A_wu = assemble_matrix(a_wu, U, W)
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)
    # A_vu = assemble_matrix(a_vu, U, V)

    AA = A_vϕ * (inv(Matrix(A_wϕ))) * A_wu
    ρV = R^2*tan(30*π/180)/2
    Mₐ = (real(AA)  )/ ρV / (ω^2) 
    Cₐ = imag(AA) / (ω^2) / ρV
    return Mₐ, Cₐ
end # function

function main_sbm(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,nΓ,W,Φ,V,U,pmid)

    a_wϕ(ϕ, w) = ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*(w*ϕ) )dΓf⁻ - ∫( im*k*ϕ*w )dΓo + ∫(w*(n(0,pmid)⋅nΓ)*((∇∇(ϕ)⋅d(0,pmid) + ∇(ϕ))⋅n(0,pmid)) - w* ∇(ϕ)⋅nΓ)dΓ # +∫( 0.1 * jump(∇(w)⋅nE⁰) * jump( ∇(ϕ)⋅nE⁰) )dE⁰
    a_vϕ(ϕ, v) = ∫( (v⋅n(0,pmid)) *(n(0,pmid)⋅nΓ) * (  im*ω*ϕ + im*ω*∇(ϕ)⋅d(0,pmid) )* (-1.0))dΓ
    a_wu(u, w) =  ∫(w * im*ω*(u⋅n(0,pmid))*(n(0,pmid)⋅nΓ) )dΓ

    A_wϕ = assemble_matrix(a_wϕ, Φ, W)
    A_wu = assemble_matrix(a_wu, U, W)
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)

    AA = A_vϕ * (inv(Matrix(A_wϕ))) * A_wu
    ρV = R^2*tan(30*π/180)/2
    # @show ρV = (12*4-∑(∫(1.0)dΩ)) #π*R^2/2 
    # stop
    Mₐ = (real(AA)  )/ ρV / (ω^2) 
    Cₐ = imag(AA) / (ω^2) / ρV
    return Mₐ, Cₐ
end # function

KRs = 0.1:0.1:2.0
mass_plots = plot(title="Added Mass vs KR", xlabel="KR [-]", ylabel="Added Mass Mₐ [-]")
damping_plots = plot(title="Added Damping vs KR", xlabel="KR [-]", ylabel="Added Damping Cₐ [-]")

# Loop over each model and compute added mass and damping
for (i, model) in enumerate(models)
    added_mass = []
    added_damping = []
    geo = geos[i]
    cutgeo, cutgeo_facets = cutting_model(model,geo)
    order = 1
    degree = 2 * order
    pmid = Point(0.0,-R)
        # agfem
    (Ω,Ω⁻act),(Γ,Γf⁻,Γi,Γo,E⁰),(nΓ,nE⁰),(dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰) = setup_interiors(model,cutgeo,cutgeo_facets,degree)
    W,Φ,V,U = setup_spaces(order, model, Ω⁻act, cutgeo)

    # cutfem
    # (Ω,Ω⁻act),(Γ,Γf⁻,Γi,Γo,E⁰),(nΓ,nE⁰),(dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰) = setup_interiors(model,cutgeo,cutgeo_facets,degree)
    # W,Φ,V,U = setup_spaces(order, model, Ω⁻act)

    # sbm
    # (Ω),(Γ,Γf⁻,Γi,Γo),(nΓ),(dΩ,dΓ,dΓf⁻,dΓi,dΓo) = setup_interiors(model,cutgeo,degree)
    # W,Φ,V,U = setup_spaces(order, model, Ω)

    for KR in KRs
        k = KR/R
        ω = √(k * g)
        Mₐ, Cₐ = main_agfem(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,nΓ,W,Φ,V,U)
        # Mₐ, Cₐ = main_cutfem(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰,nΓ,nE⁰,W,Φ,V,U)
        # Mₐ, Cₐ = main_sbm(k,ω,dΩ,dΓ,dΓf⁻,dΓi,dΓo,nΓ,W,Φ,V,U,pmid)
        push!(added_mass, Mₐ[2, 2])  # Extract relevant component for plotting
        push!(added_damping, Cₐ[2, 2])
        println(Mₐ[2, 2],",",Cₐ[2, 2])
    end
    println("=,=")
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