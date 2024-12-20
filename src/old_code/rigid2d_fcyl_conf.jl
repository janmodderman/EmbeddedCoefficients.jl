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
using GridapGmsh

# function agfem()
to = TimerOutput()
@timeit to "variables" begin
    # variables
    g = 9.81        # [kg/s²] gravitational constant
    ϕ₀(x) = 1.0 +im*0.0
end

@timeit to "domain" begin
    model1 = GmshDiscreteModel("data/meshes/cylinder_lowest.msh")
    model2 = GmshDiscreteModel("data/meshes/cylinder_low.msh")
    # model3 = GmshDiscreteModel("data/meshes/cylinder.msh")
end

function main(KR,model)
    @timeit to "embedding" begin
        Ω = Interior(model) 
        Γ = BoundaryTriangulation(model, tags=["cylinder"])
        nΓ = get_normal_vector(Γ)

        # create integration space (triangulation) & Gauss quadratures (measure)
        order = 1
        degree = 2*order
        Γf⁻ = BoundaryTriangulation(model, tags=["surface"])
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
        W = FESpace(Ω, reffeᵩ, vector_type=Vector{ComplexF64})
        Φ = TrialFESpace(W)
        E = FESpace(Ω, reffeᵩ, vector_type=Vector{ComplexF64})
        Η = TrialFESpace(E)
        V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{2,ComplexF64})
        # V = FESpace(Ω, reffeu, vector_type=Vector{ComplexF64})
        U = TrialFESpace(V)

        # final FE spaces
        X = MultiFieldFESpace([Φ, U])
        Y = MultiFieldFESpace([W, V])
    end

    R = 0.0756
    k = KR/R
    ω =  √(k*g)#  2π/1.2 +0.0*im
    nz = VectorValue(0.0,1.0)
    m = 0.96 + im*0.0
    # k = ω^2/g
    η₀ = 0.05
    ∇ₙϕd(x) = -im*ω*η₀*exp(im*k*x[2])

    # a((ϕ, u), (w, v)) =  ∫( ∇(ϕ)⋅∇(w) )dΩ - ∫(w * im*ω*(u⋅nΓ))dΓ +  ∫(w * (∇(ϕ)⋅nΓf ))dΓf⁻ + ∫(v ⋅ ( g*(u⋅nz)*nz + im*ω*ϕ - m*(ω^2)*u))dΓ
    # l((w, v)) = ∫(w*0.0)dΩ - ∫(w*∇ₙϕd)dΓf⁻
    # l₀1(w) = ∫(w*0.0)dΩ 
    # l₀2(w) = ∫(w⋅ VectorValue(0.0,0.0))dΩ 
    # # solver definition
    # ls = LUSolver()
    # op = AffineFEOperator(a,l,X,Y)

    a_wϕ(ϕ, w) = ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*(w*ϕ) )dΓf⁻ - ∫( im*k*ϕ*w )dΓi - ∫( im*k*ϕ*w )dΓo
    a_vϕ(ϕ, v) = ∫( (v⋅nΓ) * (  im*ω*ϕ )* (-1.0))dΓ
    a_wu(u, w) =  ∫(w * im*ω*(u⋅nΓ) )dΓ
    a_vu(u, v) = ∫(v ⋅ ( g*(u⋅nz)*nz + m*(ω^2)*u))dΓ

    A_wϕ = assemble_matrix(a_wϕ, Φ, W)
    A_wu = assemble_matrix(a_wu, U, W)
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)
    A_vu = assemble_matrix(a_vu, U, V)

    AA = A_vϕ*(inv(Matrix(A_wϕ)))*A_wu
    ρV = (π*R^2)
    Mₐ = real(AA)/ρV/(ω^2)
    # Mₐ = tr(real(AA)/ρV/(ω^2))
    Cₐ = imag(AA)/(ω^2)/ρV
    # Cₐ = tr(imag(AA)/(ω^2)/ρV)
    # @show Mₐ
    # @show Cₐ
    return Mₐ,Cₐ
end
KRs = 0.01:0.05:2.01
main(KRs[1],model1)
plt = plot()
lst = []
plt2 = plot()
lst2 = []
for KR in KRs
    @show vals = main(KR,model1)
    push!(lst,vals[1][2,2])
    push!(lst2,vals[2][2,2])
end
plot!(KRs, lst)
display(plt)


# for KR in KRs
#     # @show vals = main(KR,model1)
#     push!(lst2,vals[2][2,2])
# end
plot!(KRs, lst2)
display(plt2)
# main(model2)
# main(model3)


end

# end