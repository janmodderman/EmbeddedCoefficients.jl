module embedded_2d_CYLINDER
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
using DataFrames:DataFrame
using DataFrames:Matrix
using Plots

# function agfem()
to = TimerOutput()
@timeit to "variables" begin
    # variables
    g = 9.81        # [kg/s²] gravitational constant
    n₁ = 200#360        # [-] number of elements horizontally
    n₂ = 200         # [-] number of elements vertically
    println("Number of elements: ", n₁*n₂)
    D = 0.2
    B = 0.4
    h = 0.4        # [m] depth
    L₁ = 1.2
    ϕ₀(x) = 1.0 +im*0.0
end

@timeit to "bgmodel setup" begin
    pmin = Point(-L₁/2, -h)
    pmax = Point(L₁/2, 0.0)
    pmid = 0.5*(pmax + pmin) + VectorValue(0.0, h/2) - VectorValue(B/2,D/2)
    geo = quadrilateral(;x0=pmid,d1=VectorValue(B,0.0),d2=VectorValue(0.0,D))
    partition = (n₁, n₂)
    cart_model = CartesianDiscreteModel(pmin, pmax, partition)
    # bgmodel = simplexify(cart_model,positive=true)
    model = cart_model
    # assign labels to mesh surfaces and corner points
    labels_Γ = get_face_labeling(model)  # get the face labeling of model_Γ 
    add_tag_from_tags!(labels_Γ, "seabed", [1,2,5])  
    add_tag_from_tags!(labels_Γ, "inlet", [3, 7])
    add_tag_from_tags!(labels_Γ, "outlet", [4, 8]) 
    add_tag_from_tags!(labels_Γ, "surface", [6]) 
end


@timeit to "embedding" begin
    # define geometry & apply embedded boundary 
    cutgeo = cut(model, !geo)
    cutgeo_facets = cut_facets(model, !geo)
    Ω = Interior(model) 
    Ω⁻ = Interior(cutgeo, PHYSICAL)
    Ω⁻act = Interior(cutgeo, ACTIVE)
    Γ = EmbeddedBoundary(cutgeo)
    nΓ = get_normal_vector(Γ)

    # create integration space (triangulation) & Gauss quadratures (measure)
    order = 1
    degree = 2*order
    Γf⁻ = BoundaryTriangulation(cutgeo_facets, PHYSICAL, tags=["surface"])
    dΓf⁻ = Measure(Γf⁻, degree)
    nΓf = get_normal_vector(Γf⁻)
    dΩ⁻ = Measure(Ω⁻, degree)
    dΩ = Measure(Ω, degree)
    dΓ = Measure(Γ, degree)

    # boundaries
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
    # definition of FE spaces
    reffeᵩ = ReferenceFE(lagrangian, Float64, order)

    # AgFEM
    Wstd = FESpace(Ω⁻act, reffeᵩ, vector_type=Vector{ComplexF64})
    Φstd = TrialFESpace(Wstd)
    Dstd = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{2,ComplexF64})
    Rstd = TrialFESpace(Dstd)

    threshold = 1.0
    strategy = AggregateCutCellsByThreshold(threshold)
    aggregates = aggregate(strategy, cutgeo, geo, OUT)
    W = AgFEMSpace(Wstd, aggregates)
    Φ = TransientTrialFESpace(W)


    # final FE spaces
    X = MultiFieldFESpace([Φ, Rstd])
    Y = MultiFieldFESpace([W, Dstd])
end

ω =  2π/1.2 +0.0*im
nz = VectorValue(0.0,1.0)
m = 0.96 + im*0.0
k = ω^2/g
η₀ = 0.05
∇ₙϕd(x) = -im*ω*η₀*exp(im*k*x[1])

a((ϕ, u), (w, v)) =  ∫( ∇(ϕ)⋅∇(w) )dΩ⁻ - ∫(w * im*ω*(u⋅nΓ))dΓ +  ∫(w * (∇(ϕ)⋅nΓf ))dΓf⁻ + ∫(v ⋅ ( g*(u⋅nz)*nz + im*ω*ϕ - m*(ω^2)*u))dΓ
l((w, v)) = ∫(w*0.0)dΩ⁻ - ∫(w*∇ₙϕd)dΓf⁻
# solver definition
ls = LUSolver()
op = AffineFEOperator(a,l,X,Y)

# solve
# ϕhₜ = Gridap.solve(ode_solver, op, t₀, Tf, (x₀, v₀,v₀))
(ϕh,uh) = Gridap.solve(op)


folder="data/sims/hydrocoeffs"
name="AgFEM"

    writevtk(Ω⁻, folder*"/"*name*".vtu", cellfields=[ "phih_r"=>real(ϕh), "phih_i"=>imag(ϕh) ])
    writevtk(Γ, folder*"/"*name*"_u.vtu", cellfields=["uh_r"=>real(uh), "uh_i"=>imag(uh)])

show(to)
end

# end