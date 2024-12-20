module RECTANGULAR_SPECTRAL
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
# using GridapGmsh

# function agfem()
to = TimerOutput()
@timeit to "variables" begin
    # variables
    g = 9.81        # [kg/s²] gravitational constant
    n₁ = 120#360        # [-] number of elements horizontally
    n₂ = 40        # [-] number of elements vertically
    println("Number of elements: ", n₁*n₂)
    D = 0.1
    B = 2*D
    h = 40*D       # [m] depth
    # ω = √(0.1*2*g/B)
    k_max = 1.0
    λ = 2*π/k_max
    @show L₁ = 2*λ
    ϕ₀(x) = 1.0 +im*0.0
end

@timeit to "domain" begin
    # model1 = GmshDiscreteModel("data/meshes/cylinder_lowest.msh")
    # model2 = GmshDiscreteModel("data/meshes/cylinder_low.msh")
    # model3 = GmshDiscreteModel("data/meshes/cylinder.msh")
    pmin = Point(0.0, -h)
    pmax = Point(L₁, 0.0)
    pmid = Point(0.0,0.0)
    geo = quadrilateral(;x0=pmid,d1=VectorValue(B/2,0.0),d2=VectorValue(0.0,D))
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

function main(KR,model)
    
    @timeit to "embedding" begin
        # Ω = Interior(model) 
        cutgeo = cut(model, !geo)
        cutgeo_facets = cut_facets(model, !geo)
        # Ω = Interior(model) 
        Ω = Interior(cutgeo, PHYSICAL)
        Ω⁻act = Interior(cutgeo, ACTIVE)
        Γ = EmbeddedBoundary(cutgeo)
       # Γ = BoundaryTriangulation(model, tags=["cylinder"])
        nΓ = get_normal_vector(Γ)

        # create integration space (triangulation) & Gauss quadratures (measure)
        order = 1
        degree = 2*order
        Γf⁻ = BoundaryTriangulation(cutgeo_facets, PHYSICAL, tags=["surface"])#BoundaryTriangulation(model, tags=["surface"])
        dΓf⁻ = Measure(Γf⁻, degree)
        nΓf = get_normal_vector(Γf⁻)
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
        reffeu = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
        # AgFEM
        W = FESpace(Ω⁻act, reffeᵩ, vector_type=Vector{ComplexF64})
        Φ = TrialFESpace(W)

        # threshold = 1.0
        # strategy = AggregateCutCellsByThreshold(threshold)
        # aggregates = aggregate(strategy, cutgeo, geo, OUT)
        # W = AgFEMSpace(Wstd, aggregates)
        # Φ = TransientTrialFESpace(W)

        # E = FESpace(Ω, reffeᵩ, vector_type=Vector{ComplexF64})
        # Η = TrialFESpace(E)
        V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{2,ComplexF64})
        # V = FESpace(Ω, reffeu, vector_type=Vector{ComplexF64})
        U = TrialFESpace(V)

        # final FE spaces
        X = MultiFieldFESpace([Φ, U])
        Y = MultiFieldFESpace([W, V])
    end

    k = KR
    ω =  √(k*g*tanh(k*h))#  2π/1.2 +0.0*im
    @show ω2b2g = (ω^2)*B/(2*g)
    # nz = VectorValue(0.0,1.0)
    # m = 0.96 + im*0.0


    a_wϕ(ϕ, w) = ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*tanh(k*h)*(w*ϕ) )dΓf⁻ - ∫( im*k*ϕ*w )dΓo
    a_vϕ(ϕ, v) = ∫( (v⋅nΓ) * (  im*ω*ϕ )* (-1.0))dΓ
    a_wu(u, w) =  ∫(w * im*ω*(u⋅nΓ) )dΓ
    # a_vu(u, v) = ∫(v ⋅ ( g*(u⋅nz)*nz + m*(ω^2)*u))dΓ

    A_wϕ = assemble_matrix(a_wϕ, Φ, W)
    A_wu = assemble_matrix(a_wu, U, W)
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)
    # A_vu = assemble_matrix(a_vu, U, V)

    AA = A_vϕ*(inv(Matrix(A_wϕ)))*A_wu
    ρV = D*B/2
    Mₐ = real(AA)/ρV/ (ω^2)
    # Mₐ = tr(real(AA)/ρV/(ω^2))
    Cₐ = imag(AA)/ρV / (ω^2)
    # Cₐ = tr(imag(AA)/(ω^2)/ρV)
    # @show Mₐ
    # @show Cₐ

    return Mₐ,Cₐ,ω2b2g
end
KRs = 1.0:0.5:25.0
main(KRs[1],model)
plt = plot()
lst = []
plt2 = plot()
lst2 = []
lstω = []
for KR in KRs
    @show vals = main(KR,model)
    push!(lst,vals[1][2,2])
    push!(lst2,vals[2][2,2])
    push!(lstω,vals[3])
end
plot!(lstω, lst, labels=["added mass"])
display(plt)


# for KR in KRs
#     # @show vals = main(KR,model1)
#     push!(lst2,vals[2][2,2])
# end
plot!(lstω, lst2, labels=["added damping"])
display(plt2)
# main(model2)
# main(model3)


end

# end