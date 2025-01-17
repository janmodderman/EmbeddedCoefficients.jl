module case1
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
KRs = [0.1:0.1:2.5;]      # range of non-dimensional wave numbers
order = 1                   # order of elements, either 1 or 2
γg = 0.1                    # GP stabilization parameter 
h = 0.0035                  # smallest element size in background mesh
outputdir = "data/sims/"   # output directory

# case specific variables
# R = 0.1                         # [m]: radius
D = 4.0
H = 4.0 # [m]: draft
sp = 5.2
# B = 2*D
# ρV = D*B/2                    # [m]: area of a full horizontal cylinder (half domain)

Ks = KRs./(D/2)                     # [m⁻¹]: range of wave numbers
# HR_ratios = [0.000, 0.342, 0.643, 0.809, 0.906, -0.259, -0.643, -0.809, -0.906]
# HR_names = ["0000", "0342", "0643", "0809", "0906", "m0259", "m0643", "m0809", "m0906"]
HR_ratios = [0.0]#[0.0,0.342, 0.643, 0.809, 0.906]
HR_names = ["0000","0342", "0643", "0809", "0906"]
# fs = [0.01:0.025:0.8;]
# ωs = 2*π.*fs

for (i,HR) in enumerate(HR_ratios)
    name="case1"
    # load background model, and geometry; then cut geometry into model
    # model, geo = helpers.setup_domain("data/meshes/background_oc4.msh",Val(3),Val(:oc4))
    model = GmshDiscreteModel("data/meshes/background_case1.msh",orient_if_simplex=false)
    # geo = STLGeometry("data/meshes/oc4_nobraces.stl")
    geo1 = tube(D/2,H,x0=VectorValue(-sp/2,0.0,-H),v=VectorValue(0,0,1))

    geo2 = tube(D/2,H,x0=VectorValue(sp/2,0.0,-H),v=VectorValue(0,0,1))
    # model = simplexify(model1,positive=true)
    # model = CartesianDiscreteModel(Point(-100,-100,-100),Point(100,100,0),(50,50,50))
    # writevtk(model,"model")
    # labels = get_face_labeling(model)
    # add_tag_from_tags!(labels, "surface", [22])                                                            
    # add_tag_from_tags!(labels, "seabed", [21])                                                            
    # add_tag_from_tags!(labels, "walls", [23,24])                                                            
    # add_tag_from_tags!(labels, "inlet", [25])                                                            
    # add_tag_from_tags!(labels, "outlet", [26])                                                            

    # model = simplexify(model)
    cutgeo1 = cut(model,geo1)
    cutgeo2 = cut(model,geo2)
    cutgeo = cut(model, union(geo1,geo2))
    cutgeofacets = cut_facets(model,union(geo1,geo2))
    Ωᵢ = Interior(cutgeo, PHYSICAL)
    dΩᵢ = Measure(Ωᵢ,4)
    ρV = ∑(∫(1.0)dΩᵢ)
    # ρV = π*(D/2)^2*H*2
    for k in Ks
    # for ω in ωs
    # k = Ks[1]
    ω = √(k * g)
    # k = ω^2/g
    aₐ = []
    bₐ = []
    degree = 2*order

    Ω = Interior(cutgeo, PHYSICAL_OUT)
    Ω⁻act = Interior(cutgeo, ACTIVE_OUT)
    writevtk(Ω,"case1_phy")
    writevtk(Ω⁻act,"case1_act")
    Γ1 = EmbeddedBoundary(cutgeo1)
    Γ2 = EmbeddedBoundary(cutgeo2)
    writevtk(Γ1,"case1_phygam1")
    writevtk(Γ2,"case1_phygam2")

    nΓ1 = -get_normal_vector(Γ1)
    nΓ2 = -get_normal_vector(Γ2)
    Γf = BoundaryTriangulation(cutgeofacets, PHYSICAL_OUT, tags=["surface"])
    dΓf = Measure(Γf, degree)
    dΩ = Measure(Ω, degree)
    dΓ1 = Measure(Γ1, degree)
    dΓ2 = Measure(Γ2, degree)
    Γi = BoundaryTriangulation(model, tags=["inlet"])
    dΓi = Measure(Γi, degree)
    Γo = BoundaryTriangulation(model, tags=["outlet"])
    dΓo = Measure(Γo, degree)
    Γw1 = BoundaryTriangulation(model, tags=["wall1"])
    dΓw1 = Measure(Γw1, degree)
    Γw2 = BoundaryTriangulation(model, tags=["wall2"])
    dΓw2 = Measure(Γw2, degree)
    E = GhostSkeleton(cutgeo)
    dE = Measure(E,degree)
    nE = get_normal_vector(E)

    threshold = 1.0
    strategy = AggregateCutCellsByThreshold(threshold)
    # aggregates = aggregate(strategy, cutgeo)
    aggregates = aggregate(strategy, cutgeo, union(geo1,geo2), OUT)
    reffe = ReferenceFE(lagrangian, Float64, order)
    Wstd = FESpace(Ω⁻act, reffe, vector_type=Vector{ComplexF64})
    W = AgFEMSpace(Wstd,aggregates)
    Φ = TrialFESpace(W)
    V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{3, ComplexF64})
    U = TrialFESpace(V)

    V2 = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{3, ComplexF64})
    U2 = TrialFESpace(V2)

    a_wϕ = (ϕ, w) -> ∫( ∇(ϕ)⋅∇(w) )dΩ - ∫(k*(w*ϕ) )dΓf - ∫( im*k*ϕ*w )dΓo - ∫( im*k*ϕ*w )dΓw1 - ∫( im*k*ϕ*w )dΓi - ∫( im*k*ϕ*w )dΓw2
    a_vϕ = (ϕ, v) -> ∫( (v⋅nΓ1) * (  im*ω*ϕ )* (-1.0))dΓ1 
    a_wu = (u, w)-> ∫(w * im*ω*(u⋅nΓ1) )dΓ1   
    
    a_vϕ2 = (ϕ, v) -> ∫( (v⋅nΓ2) * (  im*ω*ϕ )* (-1.0))dΓ2
    a_wu2 = (u, w)-> ∫(w * im*ω*(u⋅nΓ2) )dΓ2  
    
    A_wϕ,A_wu,A_vϕ = helpers.assemble_matrices(a_wϕ,a_wu,a_vϕ,W,V,Φ,U)
    _,A_wu2,A_vϕ2 = helpers.assemble_matrices(a_wϕ,a_wu2,a_vϕ2,W,V2,Φ,U2)

    @timeit to "inverse_agfem" begin
        x = helpers.Ay(A_wϕ,A_wu,A_vϕ)
        x2 = helpers.Ay(A_wϕ,A_wu2,A_vϕ2)
    end # time
    A, B = helpers.hydro_coeffs(ω,ρV,x)
    A2, B2 = helpers.hydro_coeffs(ω,ρV,x2)

    # A = real(x)*1025/(ω^2)
    # B = imag(x)*1025/ω
    push!(aₐ, A)  
    push!(bₐ, B)

    println(k*D ," ",A," ",A2," ",B," ",B2)
    # println(k*D ," ",A," ",A2[1,1]," ",B[1,1]+B2[1,1])
    # println(k*D ," ",A2[3,3]," ",B2[3,3])
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