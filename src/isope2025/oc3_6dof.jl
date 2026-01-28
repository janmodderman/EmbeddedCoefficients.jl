module oc3
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
using GridapEmbedded.CSG
include("../helper_functions.jl")

to = TimerOutput()
# global variables      
g = 9.81                    # [kg m/s²]: gravitational constant
KRs = [0.1:0.1:2.6;]      # range of non-dimensional wave numbers
order = 1                   # order of elements, either 1 or 2
γg = 0.1                    # GP stabilization parameter 
h = 0.0035                  # smallest element size in background mesh
outputdir = "data/sims/"   # output directory

# case specific variables
# R = 0.1                         # [m]: radius
# D = 0.1
# B = 2*D
# ρV = D*B/2                    # [m]: area of a full horizontal cylinder (half domain)

Ks = KRs                     # [m⁻¹]: range of wave numbers
# HR_ratios = [0.000, 0.342, 0.643, 0.809, 0.906, -0.259, -0.643, -0.809, -0.906]
# HR_names = ["0000", "0342", "0643", "0809", "0906", "m0259", "m0643", "m0809", "m0906"]
HR_ratios = [0.0]#[0.0,0.342, 0.643, 0.809, 0.906]
HR_names = ["0000","0342", "0643", "0809", "0906"]
#fs = [0.01:0.025:0.8;]
ωs = [0.1:0.25:5.0;]
# ωs = 2*π.*fs

for (i,HR) in enumerate(HR_ratios)
    name="OC3_test"
    # load background model, and geometry; then cut geometry into model
    # model, geo = helpers.setup_domain("data/meshes/background_oc4.msh",Val(3),Val(:oc4))
    model = GmshDiscreteModel("data/meshes/backoc3test3h.msh",orient_if_simplex=false)
    # model = GmshDiscreteModel("data/meshes/background_oc3_c7.msh",orient_if_simplex=false)
    geo = STLGeometry("data/meshes/oc3.stl")

    # geo1 = tube(6.5/2,14,x0=Point(0,0,-4),v=VectorValue(0,0,1))
    # geo4 = cylinder(9.4/2,x0=Point(0,0,-120),v=VectorValue(0,0,1))

    # function _tapered_cylinder(x,R1,R0,H,x0,v)
    #     w = x - x0  # Vector from base point x0
    #     A = w ⋅ v  # Projection of w onto axis v
    #     H2 = w ⋅ w  # Squared length of w
    #     B = H2 - A*A  # Perpendicular distance squared from axis
    #     # Define a linearly varying radius R(A)
    #     R = R0 + (R1 - R0) * (A / H)  # Radius varies from R0 to R1 along height H  
    #     # Modified equation for a tapered cylinder
    #     B - R^2
    # end # function

    # function tapered_cylinder(R1,R0;x0=zero(Point{3,eltype(R)}),v=VectorValue(1,0,0),name="tapered_cylinder")
    #     H = norm(v)
    #     d = v/H
    #     function taperedcylinderfun(x)
    #         _tapered_cylinder(x,R1,R0,H,x0,d)
    #     end # function
    #     tree = Leaf((taperedcylinderfun,name,nothing))
    #     AnalyticalGeometry(tree)
    # end # function

    # geo3 = tapered_cylinder(6.5/2,9.4/2,x0=Point(0,0,-12),v=VectorValue(0,0,8))
    # geo = union(geo1,intersect(intersect(geo4,plane(x0=Point(0,0,-120),v=VectorValue(0,0,-1))),geo3))




    ###########
    # geo = sphere(30)
    # model = simplexify(model1,positive=true)
    # model = CartesianDiscreteModel(Point(-60,-60,-200),Point(60,60,0),(30,30,30))
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
    # cutfacets = cut_facets(model, geo)
    # Ωᵢ = Interior(cutgeo, PHYSICAL)
    Ωᵢ = Interior(cutgeo.cut, PHYSICAL)
    dΩᵢ = Measure(Ωᵢ,4)
    ρV = ∑(∫(1.0)dΩᵢ)
    Vt = π*(6.5/2)^2*4+π*(9.4/2)^2*108+π/3*8*((9.4/2)^2+(6.5/2)^2+(9.4/2)*(6.5/2))
    println("Theoretical volume: ", Vt, " m³; Estimated volume: ", ρV," m³")
    println("Theoretical buoyancy: ", Vt*1025*9.81, " N; Estimated buoyancy: ", ρV*1025*9.81 ," N; Reported buoyancy: ",80708100.0," N")
    # 80 708 100 N

    plt = plot(legend=:right,xlabel="ω [rad/s]",ylabel="Force [kg]")
    plt2 = plot(legend=:right,xlabel="ω [rad/s]",ylabel="Moment [kg m²]")
    plt3 = plot(legend=:right,xlabel="ω [rad/s]",ylabel="Force-Moment [kg m]")
    plt4 = plot(legend=:right,xlabel="ω [rad/s]",ylabel="Force [kg/s]")
    plt5 = plot(legend=:right,xlabel="ω [rad/s]",ylabel="Moment [kg m²/s]")
    plt6 = plot(legend=:topright,xlabel="ω [rad/s]",ylabel="Force-Moment [kg m/s]")
    aₐ = []

    aₐ1 = []
    aₐ2 = []
    aₐ3 = []
    aₐ4 = []
    aₐ5 = []
    aₐ6 = []
    aₐ15 = []
    aₐ24 = []

    bₐ1 = []
    bₐ2 = []
    bₐ3 = []
    bₐ4 = []
    bₐ5 = []
    bₐ6 = []
    bₐ15 = []
    bₐ24 = []

    bₐ = []

    # for k in Ks
    for ω in ωs
    # k = Ks[1]
    # ω = √(k * g)
    # ω = √(k * g*tanh(k*320))
    k = ω^2/g
    degree = 2*order

    Ω = Interior(cutgeo.cut, PHYSICAL_OUT)
    # Ω = Interior(cutgeo, PHYSICAL_OUT)
    Ω⁻act = Interior(cutgeo.cut, ACTIVE_OUT)
    # Ω⁻act = Interior(cutgeo, ACTIVE_OUT)
    # writevtk(Ω,"oc3_phy")
    # writevtk(Ω⁻act,"oc4_act")
    Γ = EmbeddedBoundary(cutgeo.cut)
    # Γ = EmbeddedBoundary(cutgeo)
    # writevtk(Γ,"oc3_phygam")

    nΓ = -get_normal_vector(Γ)
    # writevtk(Γ,"oc4_phygam",cellfields=["n"=>nΓ])

    Γf = BoundaryTriangulation(cutgeo.cutfacets, PHYSICAL_OUT, tags=["surface"])
    # Γf = BoundaryTriangulation(cutfacets, PHYSICAL_OUT, tags=["surface"])
    # writevtk(Γf,"oc4_fsgam")
    dΓf = Measure(Γf, degree)
    dΩ = Measure(Ω, degree)
    dΓ = Measure(Γ, degree)
    # Γi = BoundaryTriangulation(model, tags=["inlet"])
    # dΓi = Measure(Γi, degree)
    # Γo = BoundaryTriangulation(model, tags=["outlet"])
    # dΓo = Measure(Γo, degree)
    # Γw1 = BoundaryTriangulation(model, tags=["wall1"])
    # dΓw1 = Measure(Γw1, degree)
    # Γw2 = BoundaryTriangulation(model, tags=["wall2"])
    # dΓw2 = Measure(Γw2, degree)

    Γw = BoundaryTriangulation(model, tags=["walls"])
    dΓw = Measure(Γw, degree)

    E = GhostSkeleton(cutgeo)
    dE = Measure(E,degree)
    nE = get_normal_vector(E)

    threshold = 1.0
    strategy = AggregateCutCellsByThreshold(threshold)
    aggregates = aggregate(strategy, cutgeo, geo, OUT)
    reffe = ReferenceFE(lagrangian, Float64, order)
    Wstd = FESpace(Ω⁻act, reffe, vector_type=Vector{ComplexF64})
    W = AgFEMSpace(Wstd,aggregates)
    Φ = TrialFESpace(W)
    V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{6, ComplexF64})
    U = TrialFESpace(V)     
 
    ΔL = 50
    χ(x) = min(1,(70-√(x[1]^2+x[2]^2))/ΔL)
    α = x -> (1.0 + im*(1-(exp(χ(x)^3.5)-1)/(exp(1)-1)))^2

    r(x) = x
    # id_r = TensorValue{3,6}(0,0,0, 1,0,0, 0,0,0, 0,1,0, 0,0,0, 0,0,1)
    # id_t = TensorValue{3,6}(1,0,0, 0,0,0, 0,1,0 ,0,0,0, 0,0,1, 0,0,0)

    id_r = TensorValue{3,6}(0,0,0, 0,0,0, 0,0,0, 1,0,0, 0,1,0, 0,0,1)
    id_t = TensorValue{3,6}(1,0,0, 0,1,0, 0,0,1 ,0,0,0, 0,0,0, 0,0,0)

    
    nz = VectorValue(0,0,1)
    nx = VectorValue(1,0,0)
    ny = VectorValue(0,1,0)
    n = nz

    a_wϕ = (ϕ, w) -> ∫( ∇(ϕ)⋅∇(w) )dΩ - (ω^2)/g*∫( w*ϕ )dΓf - ∫( im*k*ϕ*w )dΓw #- ∫( im*k*ϕ*w )dΓo - ∫( im*k*ϕ*w )dΓw1 - ∫( im*k*ϕ*w )dΓw2  - ∫( im*k*ϕ*w )dΓi

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
    push!(aₐ1, A[1,1])  
    push!(aₐ2, A[2,2])  
    push!(aₐ3, A[3,3])  
    push!(aₐ4, A[4,4])  
    push!(aₐ5, A[5,5])  
    push!(aₐ6, A[6,6])  
    push!(aₐ15, A[1,5])  
    push!(aₐ24, A[2,4])  
    push!(bₐ, B)
    push!(bₐ1, B[1,1])  
    push!(bₐ2, B[2,2])  
    push!(bₐ3, B[3,3])  
    push!(bₐ4, B[4,4])  
    push!(bₐ5, B[5,5])  
    push!(bₐ6, B[6,6])  
    push!(bₐ15, B[1,5])  
    push!(bₐ24, B[2,4]) 

    println("Frequency: ",ω/(2π), " Hz & Radial Frequency: ", ω, " rad/s")
    println("Added mass matrix: ")
    display(A)
    println("Added damping matrix: ")
    display(B)

    helpers.write_csv(aₐ,bₐ,outputdir*"agfem/"*name*"_$order.csv";namex="A",namey="B")
end # for

data_a11=CSV.read("data/exp_pro/reference/oc3/A11.csv",DataFrame)
data_a22=CSV.read("data/exp_pro/reference/oc3/A22.csv",DataFrame)
data_a33=CSV.read("data/exp_pro/reference/oc3/A33.csv",DataFrame)
data_a44=CSV.read("data/exp_pro/reference/oc3/A44.csv",DataFrame)
data_a55=CSV.read("data/exp_pro/reference/oc3/A55.csv",DataFrame)
data_a66=CSV.read("data/exp_pro/reference/oc3/A66.csv",DataFrame)
data_a15=CSV.read("data/exp_pro/reference/oc3/A15.csv",DataFrame)
data_a24=CSV.read("data/exp_pro/reference/oc3/A24.csv",DataFrame)

data_b11=CSV.read("data/exp_pro/reference/oc3/B11.csv",DataFrame)
data_b22=CSV.read("data/exp_pro/reference/oc3/B22.csv",DataFrame)
data_b33=CSV.read("data/exp_pro/reference/oc3/B33.csv",DataFrame)
data_b44=CSV.read("data/exp_pro/reference/oc3/B44.csv",DataFrame)
data_b55=CSV.read("data/exp_pro/reference/oc3/B55.csv",DataFrame)
data_b66=CSV.read("data/exp_pro/reference/oc3/B66.csv",DataFrame)
data_b15=CSV.read("data/exp_pro/reference/oc3/B15.csv",DataFrame)
data_b24=CSV.read("data/exp_pro/reference/oc3/B24.csv",DataFrame)

plot!(plt,ωs,aₐ1,label="A11")
plot!(plt,ωs,aₐ2,label="A22")
plot!(plt,ωs,aₐ3,label="A33")
scatter!(plt1,data_a11[!,1],data_a11[!,2],label="report A11",markercolor=:black,markersize=1)
scatter!(plt1,data_a22[!,1],data_a22[!,2],label="report A22",markercolor=:black,markersize=1)
scatter!(plt1,data_a33[!,1],data_a33[!,2],label="report A33",markercolor=:black,markersize=1)

display(plt)

plot!(plt2,ωs,aₐ4,label="A44")
plot!(plt2,ωs,aₐ5,label="A55")
plot!(plt2,ωs,aₐ6,label="A66")
scatter!(plt2,data_a44[!,1],data_a44[!,2],label="report A44",markercolor=:black,markersize=1)
scatter!(plt2,data_a55[!,1],data_a55[!,2],label="report A55",markercolor=:black,markersize=1)
scatter!(plt2,data_a66[!,1],data_a66[!,2],label="report A66",markercolor=:black,markersize=1)

display(plt2)

plot!(plt3,ωs,aₐ15,label="A15")
plot!(plt3,ωs,aₐ24,label="A24")
scatter!(plt3,data_a15[!,1],data_a15[!,2],label="report A15",markercolor=:black,markersize=1)
scatter!(plt3,data_a24[!,1],data_a24[!,2],label="report A24",markercolor=:black,markersize=1)

display(plt3)

plot!(plt4,ωs,bₐ1,label="B11")
plot!(plt4,ωs,bₐ2,label="B22")
plot!(plt4,ωs,bₐ3,label="B33")
scatter!(plt4,data_b11[!,1],data_b11[!,2],label="report B11",markercolor=:black,markersize=1)
scatter!(plt4,data_b22[!,1],data_b22[!,2],label="report B22",markercolor=:black,markersize=1)
scatter!(plt4,data_b33[!,1],data_b33[!,2],label="report B33",markercolor=:black,markersize=1)
display(plt4)

plot!(plt5,ωs,bₐ4,label="B44")
plot!(plt5,ωs,bₐ5,label="B55")
plot!(plt5,ωs,bₐ6,label="B66")
scatter!(plt5,data_b44[!,1],data_b44[!,2],label="report B44",markercolor=:black,markersize=1)
scatter!(plt5,data_b55[!,1],data_b55[!,2],label="report B55",markercolor=:black,markersize=1)
scatter!(plt5,data_b66[!,1],data_b66[!,2],label="report B66",markercolor=:black,markersize=1)

display(plt5)

plot!(plt6,ωs,bₐ15,label="B15")
plot!(plt6,ωs,bₐ24,label="B24")
scatter!(plt6,data_b15[!,1],data_b15[!,2],label="report B15",markercolor=:black,markersize=1)
scatter!(plt6,data_b24[!,1],data_b24[!,2],label="report B24",markercolor=:black,markersize=1)

display(plt6)
end # for

# for (i,HR) in enumerate([0.0])
#     helpers.plotter_case41(HR_names[i],HR,order)
# end # for
end # module