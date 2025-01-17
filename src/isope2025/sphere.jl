module hemisphere
using Gridap, Gridap.Arrays, Gridap.FESpaces, Gridap.Fields, Gridap.Geometry, Gridap.TensorValues
using GridapEmbedded, GridapEmbedded.Interfaces, GridapEmbedded.LevelSetCutters
using GridapGmsh
using WriteVTK
using TimerOutputs
using DrWatson
using Plots
using LinearAlgebra
using CSV, DataFrames
using STLCutters
include("../helper_functions.jl")

to = TimerOutput()
# global variables      
g = 9.81                    # [kg m/s²]: gravitational constant
KRs = [0.1:0.1:4;]      # range of non-dimensional wave numbers
order = 1                   # order of elements, either 1 or 2
γg = 0.1                    # GP stabilization parameter 
h = 0.0035                  # smallest element size in background mesh
outputdir = "data/sims/"   # output directory

# case specific variables
# R = 0.1                         # [m]: radius
# D = 0.1
# B = 2*D
# ρV = D*B/2                    # [m]: area of a full horizontal cylinder (half domain)

Ks = KRs./30                     # [m⁻¹]: range of wave numbers
# HR_ratios = [0.000, 0.342, 0.643, 0.809, 0.906, -0.259, -0.643, -0.809, -0.906]
# HR_names = ["0000", "0342", "0643", "0809", "0906", "m0259", "m0643", "m0809", "m0906"]
HR_ratios = [0.0]#[0.0,0.342, 0.643, 0.809, 0.906]
HR_names = ["0000","0342", "0643", "0809", "0906"]
fs = [0.01:0.025:0.8;]
ωs = 2*π.*fs

for (i,HR) in enumerate(HR_ratios)
    name="sphere"
    # load background model, and geometry; then cut geometry into model
    # model, geo = helpers.setup_domain("data/meshes/background_oc4.msh",Val(3),Val(:oc4))
    model = GmshDiscreteModel("data/meshes/background_oc4_14.msh",orient_if_simplex=false)
    # geo = STLGeometry("data/meshes/oc4.stl")
    geo = sphere(30)
    # model = simplexify(model1,positive=true)
    # model = CartesianDiscreteModel(Point(0,-60,-200),Point(60,60,0),(20,20,100))
    # # writevtk(model,"model")
    # # model = simplexify(model,positive=true)
    # labels = get_face_labeling(model)
    # add_tag_from_tags!(labels, "surface", [22])                                                            
    # add_tag_from_tags!(labels, "seabed", [21])                                                            
    # add_tag_from_tags!(labels, "wall1", [23])                                                            
    # add_tag_from_tags!(labels, "wall2", [24])                                                            
    # add_tag_from_tags!(labels, "inlet", [25])                                                            
    # add_tag_from_tags!(labels, "outlet", [26])                                                            

    # model = simplexify(model)
    cutgeo = cut(model, geo)
    cutfacets = cut_facets(model,geo)
    Ωᵢ = Interior(cutgeo, PHYSICAL)
    dΩᵢ = Measure(Ωᵢ,4)
    ρV = ∑(∫(1.0)dΩᵢ)
    @show ρV*g*1025
    for k in Ks
    # for ω in ωs
    # k = Ks[1]
    ω = √(k * g*tanh(k*200))
    # k = ω^2/g
    aₐ = []
    bₐ = []
    degree = 2*order

    Ω = Interior(cutgeo, PHYSICAL_OUT)
    Ω⁻act = Interior(cutgeo, ACTIVE_OUT)
    # writevtk(Ω,"oc4_phy")
    # writevtk(Ω⁻act,"oc4_act")
    Γ = EmbeddedBoundary(cutgeo)
    # writevtk(Γ,"oc4_phygam")

    nΓ = -get_normal_vector(Γ)
    # writevtk(Γ,"oc4_phygam",cellfields=["n"=>nΓ])

    Γf = BoundaryTriangulation(cutfacets, PHYSICAL_OUT, tags=["surface"])
    # writevtk(Γf,"oc4_fsgam")
    dΓf = Measure(Γf, degree)
    dΩ = Measure(Ω, degree)
    dΓ = Measure(Γ, degree)
    Γw = BoundaryTriangulation(model, tags=["walls"])
    dΓw = Measure(Γw, degree)
    E = GhostSkeleton(cutgeo)
    dE = Measure(E,degree)
    nE = get_normal_vector(E)

    threshold = 1.0
    strategy = AggregateCutCellsByThreshold(threshold)
    # aggregates = aggregate(strategy, cutgeo)
    aggregates = aggregate(strategy, cutgeo, geo, OUT)
    reffe = ReferenceFE(lagrangian, Float64, order)
    Wstd = FESpace(Ω⁻act, reffe, vector_type=Vector{ComplexF64})
    # W = FESpace(Ω⁻act, reffe, vector_type=Vector{ComplexF64})
    W = AgFEMSpace(Wstd,aggregates)
    Φ = TrialFESpace(W)
    V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{3, ComplexF64})
    U = TrialFESpace(V)
    # r(x) = x
    # r = x->r(x)
    a_wϕ = (ϕ, w) -> ∫( ∇(ϕ)⋅∇(w) )dΩ - (ω^2)/g*∫( w*ϕ )dΓf  - ∫( im*k*ϕ*w )dΓw
    a_vϕ = (ϕ, v) -> ∫( ((id_t⋅v)⋅nΓ) * (  im*ω*ϕ )* (-1.0))dΓ + ∫( ((id_r⋅v)⋅(r×nΓ)) * (  im*ω*ϕ )* (-1.0))dΓ
    a_wu = (u, w)-> ∫(w * im*ω*((id_t⋅u)⋅nΓ) )dΓ + ∫(w * im*ω*((id_r⋅u)⋅(r×nΓ)) )dΓ   
    
    A_wϕ,A_wu,A_vϕ = helpers.assemble_matrices(a_wϕ,a_wu,a_vϕ,W,V,Φ,U)

    @timeit to "inverse_agfem" begin
        x = helpers.Ay(A_wϕ,A_wu,A_vϕ)
    end # time
    A1, B1 = helpers.hydro_coeffs(ω,ρV,x)
    A = real(x)/(ω^2)*1025  # [kg]: added mass
    B = imag(x)/ω*1025      # [kg/s]: added damping
    push!(aₐ, A)  
    push!(bₐ, B)

    println("Frequency: ",ω/(2π), " Hz & Radial Frequency: ", ω, " rad/s")
    println("Added mass matrix: ")
    display(A)
    println("Added damping matrix: ")
    display(B)

    # println("Added mass coefficient matrix: ")
    # display(A1)
    # println("Added damping coefficient matrix: ")
    # display(B1)
    # println(A[3,3])
    # println(B[3,3])

    # run case for agfem, cutfem or sbm
    # (aₐ,bₐ) = helpers.run_agfem(Ks, ρV,g, order, model, cutgeo, cutgeo_facets,to)
    # (aₑ,bₑ) = helpers.run_cutfem(Ks, ρV,g, order, model, cutgeo, cutgeo_facets, γg, h,to)
    # (Ω,_,_,_,_),(nΓ,),(dΩ,dΓ,dΓf,_,dΓo) = helpers.setup_interiors(model,cutgeo,2)
    # ρV2 = (12*4-∑(∫(1.0)dΩ))
    # (aₛ,bₛ) = helpers.run_sbm(Ks, ρV2,g, order, model, cutgeo, helpers.n(pcor,Val(:rectangle)), helpers.d(pcor,Val(:rectangle)),to)

    helpers.write_csv(aₐ,bₐ,outputdir*"agfem/"*name*"_$order.csv";namex="A",namey="B")
    # helpers.write_csv(aₑ,bₑ,outputdir*"cutfem/"*name*"_$order.csv";namex="A",namey="B")
    # helpers.write_csv(aₛ,bₛ,outputdir*"sbm/"*name*"_$order.csv";namex="A",namey="B")

end # for
end # for

# for (i,HR) in enumerate([0.0])
#     helpers.plotter_case41(HR_names[i],HR,order)
# end # for
end # module