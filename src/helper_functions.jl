 module helpers 
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

 #==============================================================================================================================#

# DOMAIN SETUP 
# horizontal cylinder (2D)
function setup_domain(R::Float64,pmid::VectorValue,mesh::String,::Val{2},::Val{:cylinder})
    model = GmshDiscreteModel(mesh)
    geo = disk(R, x0=pmid)
    model,geo
end # function

function setup_domain(R::Float64,pmid::VectorValue,L₁::Float64,d::Float64,partition::Tuple,::Val{2},::Val{:cylinder})
    pmin = Point(0.0, -d)
    pmax = Point(L₁, 0.0)
    model = CartesianDiscreteModel(pmin, pmax, partition)
    geo = disk(R, x0=pmid)
    labels_Γ = get_face_labeling(model)  # get the face labeling of model_Γ 
    add_tag_from_tags!(labels_Γ, "seabed", [1,2,5])  
    add_tag_from_tags!(labels_Γ, "inlet", [3, 7])
    add_tag_from_tags!(labels_Γ, "outlet", [4, 8]) 
    add_tag_from_tags!(labels_Γ, "surface", [6]) 
    model,geo
end # function

# rectangle (2D)
function setup_domain(D::Float64,B::Float64,pmin::VectorValue,mesh::String,::Val{2},::Val{:rectangle})
    model = GmshDiscreteModel(mesh)
    geo = quadrilateral(;x0=pmin-VectorValue(B/2,D),d1=VectorValue(B,0.0),d2=VectorValue(0.0,2*D))
    model,geo
end # function

function setup_domain(D::Float64,B::Float64,pmid::VectorValue,L₁::Float64,d::Float64,partition::Tuple,::Val{2},::Val{:rectangle})
    pmin = Point(0.0, -d)
    pmax = Point(L₁, 0.0)
    model = CartesianDiscreteModel(pmin, pmax, partition)
    geo = quadrilateral(;x0=pmid-VectorValue(B/2,D),d1=VectorValue(B,0.0),d2=VectorValue(0.0,2*D))
    labels_Γ = get_face_labeling(model)  # get the face labeling of model_Γ 
    add_tag_from_tags!(labels_Γ, "seabed", [1,2,5])  
    add_tag_from_tags!(labels_Γ, "inlet", [3, 7])
    add_tag_from_tags!(labels_Γ, "outlet", [4, 8]) 
    add_tag_from_tags!(labels_Γ, "surface", [6]) 
    model,geo
end # function

# triangle (2D)
function setup_domain(D::Float64,α::Float64,pmid::VectorValue,mesh::String,::Val{2},::Val{:triangle})
    φ=α*π/180
    geo = quadrilateral(;x0=pmid,d1=VectorValue(D*tan(φ),D),d2=VectorValue(0.0,D))
    model = GmshDiscreteModel(mesh)
    model,geo
end # function

# vertical cylinder (3D)
function setup_domain(R::Float64,D::Float64,pmid::VectorValue,mesh::String,::Val{3},::Val{:cylinder})
    model = GmshDiscreteModel(mesh)
    geo = cylinder(R;x0=pmid-VectorValue(0.0,0.0,-D),v=VectorValue(0,0,1))
    model,geo
end # function

# function setup_domain(R::Float64,D::Float64,pmid::VectorValue,L₁::Float64,d::Float64,partition::Tuple,::Val{3},::Val{:cylinder})
#     pmin = Point(-L₁/2, -L₁/2, -d)
#     pmax = Point(L₁/2, L₁/2, 0.0)
#     model = CartesianDiscreteModel(pmin, pmax, partition)
#     geo = cylinder(R;x0=pmid-VectorValue(0.0,0.0,-D),v=VectorValue(0,0,1))
#     labels_Γ = get_face_labeling(model)  # get the face labeling of model_Γ 
#     add_tag_from_tags!(labels_Γ, "seabed", [21])  
#     add_tag_from_tags!(labels_Γ, "inlet", [3, 7])
#     add_tag_from_tags!(labels_Γ, "outlet", [4, 8]) 
#     add_tag_from_tags!(labels_Γ, "surface", [22]) 
#     add_tag_from_tags!(labels_Γ, "walls", [22]) 
#     model,geo
# end # function



function setup_domain(mesh::String,::Val{3},::Val{:oc4})
    model = GmshDiscreteModel(mesh)
    geo = STLGeometry("data/meshes/oc4_nobraces.stl")
    model,geo
end # function

#==============================================================================================================================#

# CUTTING MODEL
function cutting_model(model::DiscreteModel,geo)
    cut(model, geo), cut_facets(model, geo)
end # function

function cutting_model(model::DiscreteModel,geo::STLGeometry)
    cutgeo = cut(model, geo)
    cutgeo, cutgeo.cutfacets
end # function

#==============================================================================================================================#

# FE SPACES
# cell aggregation
function setup_spaces(order::Int64, model::DiscreteModel, Ω::Triangulation, cutgeo::EmbeddedDiscretization, dim::Int64)
    threshold = 1.0
    strategy = AggregateCutCellsByThreshold(threshold)
    aggregates = aggregate(strategy, cutgeo, cutgeo.geo, OUT)
    reffe = ReferenceFE(lagrangian, Float64, order)
    Wstd = FESpace(Ω, reffe, vector_type=Vector{ComplexF64})
    W = AgFEMSpace(Wstd,aggregates)
    Φ = TrialFESpace(W)
    V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{dim, ComplexF64})
    U = TrialFESpace(V)
    W,Φ,V,U
end # function

# active/surrogate domain
function setup_spaces(order::Int64, model::DiscreteModel, Ω::Triangulation, dim::Int64)
    reffe = ReferenceFE(lagrangian, Float64, order)
    W = FESpace(Ω, reffe, vector_type=Vector{ComplexF64})
    Φ = TrialFESpace(W)
    V = ConstantFESpace(model; vector_type=Vector{ComplexF64}, field_type=VectorValue{dim, ComplexF64})
    U = TrialFESpace(V)
    W,Φ,V,U
end # function

#==============================================================================================================================#

# SETUP TRIANGULATIONS AND MEASURES
# Cut triangulations + quads
function setup_interiors(model::DiscreteModel,cutgeo::EmbeddedDiscretization,cutgeo_facets::EmbeddedFacetDiscretization,degree::Int64)
    Ω = Interior(cutgeo, PHYSICAL_OUT)
    Ω⁻act = Interior(cutgeo, ACTIVE_OUT)
    Γ = EmbeddedBoundary(cutgeo)
    nΓ = -get_normal_vector(Γ)
    Γf = BoundaryTriangulation(cutgeo_facets, tags=["surface"])
    dΓf = Measure(Γf, degree)
    dΩ = Measure(Ω, degree)
    dΓ = Measure(Γ, degree)
    Γi = BoundaryTriangulation(model, tags=["inlet"])
    dΓi = Measure(Γi, degree)
    Γo = BoundaryTriangulation(model, tags=["outlet"])
    dΓo = Measure(Γo, degree)
    E = GhostSkeleton(cutgeo)
    dE = Measure(E,degree)
    nE = get_normal_vector(E)
    (Ω,Ω⁻act,Γ,Γf,Γi,Γo,E),(nΓ,nE),(dΩ,dΓ,dΓf,dΓi,dΓo,dE)
end # function

# Surrogate triangulations + quads
function setup_interiors(model::DiscreteModel,cutgeo::EmbeddedDiscretization,degree::Int64)
    Ω = Interior(cutgeo, OUT)
    Γ = Interface(Interior(cutgeo,ACTIVE),Ω).⁻
    dΩ = Measure(Ω, degree)
    dΓ = Measure(Γ, degree)   
    nΓ = get_normal_vector(Γ)
    Γf = BoundaryTriangulation(Ω, tags=["surface"])
    dΓf = Measure(Γf, degree)
    Γi = BoundaryTriangulation(model, tags=["inlet"])
    dΓi = Measure(Γi, degree)
    Γo = BoundaryTriangulation(model, tags=["outlet"])
    dΓo = Measure(Γo, degree)
    (Ω,Γ,Γf,Γi,Γo),(nΓ,),(dΩ,dΓ,dΓf,dΓi,dΓo)
end # function

#==============================================================================================================================#

# SETUP MATRIX CONTRIBUTIONS
# AgFEM weak form := conformal weak form
function weak_form(k::Float64,ω::Float64,nΓ::CellField,dΩ::Measure,dΓ::Measure,dΓf::Measure,dΓo::Measure)
    a_wϕ = (ϕ, w) -> ∫( ∇(ϕ)⋅∇(w) )dΩ -  ∫(k*(w*ϕ) )dΓf - ∫( im*k*ϕ*w )dΓo 
    a_vϕ = (ϕ, v) -> ∫( (v⋅nΓ) * (  im*ω*ϕ )* (-1.0))dΓ
    a_wu = (u, w)-> ∫(w * im*ω*(u⋅nΓ) )dΓ
    a_wϕ,a_vϕ,a_wu
end # function

# CutFEM weak form := conformal weak form + GP terms
# GP first order
function weak_form(nE,dE,γg,h,::Val{1})
    (w,ϕ) -> ∫( (γg*(h^3)) * (jump(nE⋅∇(w)) ⊙ jump(nE⋅∇(ϕ))) )dE
end # function
# GP second order
function weak_form(nE,dE,γg,h,::Val{1})
    (w,ϕ) -> ∫( (γg*(h^3)) * (jump(nE⋅∇(w)) ⊙ jump(nE⋅∇(ϕ))) + (γg*(h^5)) * jump(nE⋅∇∇(w)) ⊙ jump(nE⋅∇∇(ϕ)) )dE
end # function

# SBM weak form := conformal weak form + shifted Neumann + shifted hydrodynamic BC
function weak_form(ω::Float64,nΓ::CellField,dΓ::Measure,n::Function,d::Function,dcf::FEFunction)
    a_wϕ = (ϕ, w) -> ∫(w*(n⋅nΓ)*((∇∇(ϕ)⋅d + ∇(ϕ))⋅n) - w* ∇(ϕ)⋅nΓ)dΓ
    a_vϕ = (ϕ, v) -> ∫( (im*ω* (-1.0))*( ϕ + (∇(ϕ)⋅d)) * (v⋅n) *(nΓ⋅n) * J(dcf))dΓ
    a_wu = (u, w) -> ∫( im*ω*w*(n⋅nΓ)*(u⋅n) )dΓ
    a_wϕ,a_vϕ,a_wu
end # function

#==============================================================================================================================#

# Jacobian J(u) for the distance function d for SBM
Ft(∇u) = TensorValue(1.0,0.0,0.0,1.0) + ∇u
J(u) = (Ft(∇(u)))⊙(Ft(∇(u)))
# J(u) = (Ft(∇(u)))⊙(Ft(∇(u)))

# pmid1(pmid,x) = pmid - x
# pmid1(pmid) = x -> pmid1(pmid,x)


# TRUE DISTANCE & NORMAL FUNCTIONS (SBM)
# horizontal cylinder (2D)
n(x,pmid,::Val{:cylinder}) = (pmid-x)/√((pmid-x)⋅(pmid-x))
# n(x,pmid) = (pmid.-x)/√((pmid.-x)⋅(pmid.-x))
# n(x,pmid) = VectorValue(1.0,0.0)#(pmid.-x)/√((pmid.-x)⋅(pmid.-x))
d(x,pmid,R,::Val{:cylinder}) =  (√((pmid-x)⋅(pmid-x))-R)*n(x,pmid)
# d(x,pmid,R) =  (√((pmid.-x)⋅(pmid.-x))-R)*n(x,pmid)
# d(x,pmid,R) = pmid(x).-x#(√((pmid.-x)⋅(pmid.-x))-R)*n(x,pmid)
n(pmid,::Val{:cylinder}) = x -> n(x,pmid,Val(:cylinder))
d(pmid,R,::Val{:cylinder}) = x -> d(x,pmid,R,Val(:cylinder))

function d(x,pcor,::Val{:rectangle})
    dx = maximum([pcor[1]-x[1], 0.0, x[1]-pcor[2]])
    dy = maximum([pcor[3]-x[2], 0.0, x[2]-pcor[4]])
    VectorValue(-dx,dy)
end # function

function n(x,pcor,::Val{:rectangle})
    dx = maximum([pcor[1]-x[1], 0.0, x[1]-pcor[2]])
    dy = maximum([pcor[3]-x[2], 0.0, x[2]-pcor[4]])
    dist = √(dx^2+dy^2)
    VectorValue(-dx/dist,dy/dist)
end # function

d(pcor,::Val{:rectangle}) = x -> d(x,pcor,Val(:rectangle))
n(pcor,::Val{:rectangle}) = x -> n(x,pcor,Val(:rectangle))

# ∇d(x,pmid,R) = 
#==============================================================================================================================#


# ASSEMBLE MATRICES
function assemble_matrices(a_wϕ::Function,a_wu::Function,a_vϕ::Function,W::FESpace,V::FESpace,Φ::FESpace,U::FESpace)
    assemble_matrix(a_wϕ, Φ, W), assemble_matrix(a_wu, U, W), assemble_matrix(a_vϕ, Φ, V)
end # function

#==============================================================================================================================#

# MATRICES PRODUCT
# very slow DO NOT USE
# function AAA(A_wϕ,A_wu,A_vϕ)
#     A_vϕ * (inv(Matrix(A_wϕ))) * A_wu
# end # function

# fast
function Ay(A_wϕ,A_wu,A_vϕ)
    sz  = size(A_wu)
    y = zeros(ComplexF64,sz[2],sz[2])
    ss = symbolic_setup(LUSolver(),A_wϕ)
    ns = numerical_setup(ss,A_wϕ)
    u = similar(Vector{ComplexF64}(A_wu[:,1]))
    for i in 1:sz[2]
        u = Gridap.Algebra.solve!(u,ns,Vector{ComplexF64}(A_wu[:,i]))
        y[i,:] = A_vϕ*u 
    end
    y
end # function

# function Ay(A_wϕ,A_wu,A_vϕ)
#     sz  = size(A_wu)
#     y = zeros(ComplexF64,sz[2],sz[2])
#     # ss = symbolic_setup(LUSolver(),A_wϕ)
#     # ns = numerical_setup(ss,A_wϕ)
#     u = similar(Vector{ComplexF64}(A_wu[:,1]))
#     for i in 1:sz[2]
#         u = Gridap.Algebra.solve!(LUSolver(),A_wϕ,Vector{ComplexF64}(A_wu[:,i]))
#         y[i,:] = A_vϕ*u 
#     end
#     y
# end # function

#==============================================================================================================================#

# HYDRODYNAMIC COEFFICIENTS
# obtain non-dimensionalized hydrodynamic coefficients
function hydro_coeffs(ω,ρV,AA)
    Mₐ = real(AA)  / ρV / (ω^2) 
    Cₐ = imag(AA) / (ω^2) / ρV
    Mₐ, Cₐ
end # function

#==============================================================================================================================#

# RUN CASES
function run_agfem(Ks, ρV,g, order, model, cutgeo, cutgeo_facets,to)
    added_mass = []
    added_damping = []
    degree = 2*order

    (_,Ω⁻act,_,_,_,_,_),(nΓ,_),(dΩ,dΓ,dΓf,_,dΓo,_) = setup_interiors(model,cutgeo,cutgeo_facets,degree)
    W,Φ,V,U = setup_spaces(order, model, Ω⁻act, cutgeo,num_dims(model))
    for k in Ks
        ω = √(k * g)
        a_wϕ,a_vϕ,a_wu = weak_form(k,ω,nΓ,dΩ,dΓ,dΓf,dΓo)
        A_wϕ,A_wu,A_vϕ = assemble_matrices(a_wϕ,a_wu,a_vϕ,W,V,Φ,U)
        @timeit to "inverse_agfem" begin
            x = Ay(A_wϕ,A_wu,A_vϕ)
        end # time
        A, B = hydro_coeffs(ω,ρV,x)
        push!(added_mass, A)  
        push!(added_damping, B)
        show(to)
    end # for
    (added_mass,added_damping)
end # function

function run_cutfem(Ks, ρV,g, order, model, cutgeo, cutgeo_facets, γg, h,to;GPflag=true)
    added_mass = []
    added_damping = []
    degree = 2*order

    (_,Ω⁻act,_,_,_,_,_),(nΓ,nE),(dΩ,dΓ,dΓf,_,dΓo,dE) = setup_interiors(model,cutgeo,cutgeo_facets,degree)
    W,Φ,V,U = setup_spaces(order, model, Ω⁻act,num_dims(model))
    for k in Ks
        ω = √(k * g)
        a_wϕ,a_vϕ,a_wu = weak_form(k,ω,nΓ,dΩ,dΓ,dΓf,dΓo)            # conformal weak form
        A_wϕ,A_wu,A_vϕ = assemble_matrices(a_wϕ,a_wu,a_vϕ,W,V,Φ,U)  # assemble matrices
        if GPflag
            a_wϕₑ = weak_form(nE,dE,γg,h,Val(order))                # GP terms
            A_wϕₑ = assemble_matrix(a_wϕₑ, Φ, W)                    # assemble GP terms matrix
            A_wϕ=A_wϕ+A_wϕₑ                                         # sum matrix contributions
        end # if
        @timeit to "inverse_cutfem" begin
            x = Ay(A_wϕ,A_wu,A_vϕ)
        end # time
        A, B = hydro_coeffs(ω,ρV,x)
        push!(added_mass, A)  
        push!(added_damping, B)
        show(to)
    end # for
    (added_mass,added_damping)
end # function

function run_sbm(Ks, ρV,g, order, model, cutgeo, n, d, to)
    added_mass = []
    added_damping = []
    degree = 2*order

    (Ω,_,_,_,_),(nΓ,),(dΩ,dΓ,dΓf,_,dΓo) = setup_interiors(model,cutgeo,degree)
    W,Φ,V,U = setup_spaces(order, model, Ω, num_dims(model))

    # required now to obtain gradient of d, TODO: find more better and elegant solution
    Vd= FESpace(Ω,ReferenceFE(lagrangian,VectorValue{num_dims(model),Float64},3)) # current implementation requires higher order to correctly get the gradient of the distance function
    dcf = interpolate_everywhere(CellField(d,Ω),Vd)
    # @show ρV2 = 1.04*(12*4-∑(∫(1.0)dΩ))
    # @show ρV3 = (∑(∫(1.0)dΓ)/2)^2
    # @show ρV4 = (∑(∫(J(dcf))dΓ)/2)^2
    for k in Ks
        ω = √(k * g)
        a_wϕ,_,_ = weak_form(k,ω,nΓ,dΩ,dΓ,dΓf,dΓo)                          # conformal weak form (only a_wϕ)
        a_wϕₛ,a_vϕ,a_wu = weak_form(ω,nΓ,dΓ,n,d,dcf)                            # shifted contributions (on a_wϕ, a_vϕ and a_uw)
        A_wϕ = assemble_matrix(a_wϕ, Φ, W)                                  # assemble conformal matrix contributions
        # A_wϕₛ,A_wu,A_vϕ = assemble_matrices(a_wϕₛ,a_wu,a_vϕ,W,V,Φ,U)        # assemble shifted matrix contributions
        A_wϕₛ = assemble_matrix(a_wϕₛ, Φ, W) 
        A_wu = assemble_matrix(a_wu, U, W) 
        A_vϕ = assemble_matrix(a_vϕ, Φ, V)
        A_wϕ=A_wϕ+A_wϕₛ                                                     # sum matrix contributions
        @timeit to "inverse_sbm" begin
            x = Ay(A_wϕ,A_wu,A_vϕ)                                              # solve for inverse
        end # time
        A, B = hydro_coeffs(ω,ρV,x)
        push!(added_mass, A)  
        push!(added_damping, B)
        # @show A[2,2]
        show(to)
    end # for
    (added_mass,added_damping)
end # function

# function run_case()

# end # function


# TODO: UNFINISHED
function plotter(x::Vector,results::String;name="")
    data1 = CSV.read("data/exp_raw/cyl/cutfemHR"*name*".csv",DataFrame)
    data2 = CSV.read("data/exp_raw/cyl/agfemHR"*name*".csv",DataFrame)
    data3 = CSV.read("data/exp_raw/cyl/sbmHR"*name*".csv",DataFrame)
    data4 = CSV.read("data/exp_raw/cyl/sbmp2HR"*name*".csv",DataFrame)
    refdatam = CSV.read("data/exp_pro/ref_cyl_me_HR"*name*".csv",DataFrame)
    refdatad = CSV.read("data/exp_pro/ref_cyl_de_HR"*name*".csv",DataFrame)

    plot!(mass_plots, KRs, data1[!,1], label="CutFEM",linecolor="#332288")#,markershape=:circle)
    plot!(mass_plots, KRs, data2[!,1], label="AgFEM",linecolor="#CC6677")#,markershape=:cross)
    plot!(mass_plots, KRs, fp1[i].*data3[!,1], label="SBM p=1",linecolor="#117733")#,markershape=:diamond)
    plot!(mass_plots, KRs, fp2[i].*data4[!,1], label="SBM p=2",linecolor="#117733",linestyle=:dash,markershape=:utriangle,markercolor="#117733")#,markershape=:diamond)

    plot!(mass_plots,refdatam[!,1],refdatam[!,2],label="h/R=0."*name[2:end],linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)
    
    plot!(damping_plots, KRs, data1[!,2], label="CutFEM",linecolor="#332288")#,markershape=:circle)
    plot!(damping_plots, KRs, data2[!,2], label="AgFEM",linecolor="#CC6677")#,markershape=:cross)
    plot!(damping_plots, KRs, fp1[i].*data3[!,2], label="SBM p=1",linecolor="#117733")#,markershape=:diamond)
    plot!(damping_plots, KRs, fp2[i].*data4[!,2], label="SBM p=2",linecolor="#117733",linestyle=:dash,markershape=:utriangle,markercolor="#117733")#,markershape=:diamond)
    
    plot!(damping_plots,refdatad[!,1],refdatad[!,2],label="h/R=0."*name[2:end],linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)
end # function

# TODO: ADD COLOR VAR
function plotter(x::Vector,y1::Vector,y2::Vector;name="")
    mass_plots = plot(xlabel="k̄ [-]", ylabel="Ā₃₃ [-]")
    damping_plots = plot(xlabel="k̄ [-]", ylabel="B̄₃₃ [-]")
    plot!(mass_plots, x, y1, label=name,linecolor="#CC6677")#,markershape=:cross)
    plot!(damping_plots, x, y2, label=name,linecolor="#CC6677")#,markershape=:cross)
    display(mass_plots)
    display(damping_plots)
end # function

function write_csv(x,y,dirname;namex="x",namey="y")
    if !isfile(dirname)
        touch(dirname)
    end # if
    flattened = [(vec(i)..., vec(j)...) for (i, j) in zip(x, y)]
    CSV.write(dirname,(DataFrame(flattened)))
end # function

# TESTING
function plotter_case42(name::String,val::Float64,order::Int64)
KRs = [0.1:0.025:2.0;]      # range of non-dimensional wave numbers

mass_plots = plot(xlabel="k̄ [-]", ylabel="Ā₃₃ [-]")
damping_plots = plot(xlabel="k̄ [-]", ylabel="B̄₃₃ [-]")
data1 = CSV.read("data/sims/cutfem/cylHR$(name)_$order.csv",DataFrame)
data2 = CSV.read("data/sims/agfem/cylHR$(name)_$order.csv",DataFrame)
data3 = CSV.read("data/sims/sbm/cylHR$(name)_$order.csv",DataFrame)
# 
iA = 4
iB = 8
msize = 4
plot!(mass_plots, KRs, data1[!,iA], label="CutFEM p = $order",linecolor="#332288")#,markershape=:circle)
plot!(mass_plots, KRs, data2[!,iA], label="AgFEM p = $order",linecolor="#CC6677")#,markershape=:cross)
plot!(mass_plots, KRs, data3[!,iA], label="SBM p = $order",linecolor="#117733")#,markershape=:diamond)
refdatam = CSV.read("data/exp_pro/ref_cyl_me_HR$name.csv",DataFrame)
refdatad = CSV.read("data/exp_pro/ref_cyl_de_HR$name.csv",DataFrame)

plot!(mass_plots,refdatam[!,1],refdatam[!,2],label="h/R=$val",linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)

plot!(damping_plots, KRs, data1[!,iB], label="CutFEM p = $order",linecolor="#332288")#,markershape=:circle)
plot!(damping_plots, KRs, data2[!,iB], label="AgFEM p = $order",linecolor="#CC6677")#,markershape=:cross)
plot!(damping_plots, KRs, data3[!,iB], label="SBM p = $order",linecolor="#117733")#,markershape=:diamond)
plot!(damping_plots,refdatad[!,1],refdatad[!,2],label="h/R=$val",linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)
display(mass_plots)
display(damping_plots)
savefig(mass_plots,"plots/case42/A_$(name)_$order.png")
savefig(damping_plots,"plots/case42/B_$(name)_$order.png")
end # function

function plotter_case41(name::String,val::Float64,order::Int64)
    KRs = [0.1:0.025:2.5;]      # range of non-dimensional wave numbers
    
    mass_plots = plot(xlabel="k̄ [-]", ylabel="Ā₃₃ [-]")
    damping_plots = plot(xlabel="k̄ [-]", ylabel="B̄₃₃ [-]")
    data1 = CSV.read("data/sims/cutfem/rectHR$(name)_$order.csv",DataFrame)
    data2 = CSV.read("data/sims/agfem/rectHR$(name)_$order.csv",DataFrame)
    data3 = CSV.read("data/sims/sbm/rectHR$(name)_$order.csv",DataFrame)
    # 
    iA = 4
    iB = 8
    msize = 4
 
    # reference data
    # experimental
    data_me = CSV.read("data/exp_pro/ref_rect_me.csv",DataFrame)
    data_de = CSV.read("data/exp_pro/ref_rect_de.csv",DataFrame)

    # HPC
    data_mhpc = CSV.read("data/exp_pro/ref_rect_mhpc.csv",DataFrame)
    data_dhpc = CSV.read("data/exp_pro/ref_rect_dhpc.csv",DataFrame)

    # linear FEM
    data_mfem1 = CSV.read("data/exp_pro/ref_rect_mfem1.csv",DataFrame)
    data_dfem1 = CSV.read("data/exp_pro/ref_rect_dfem1.csv",DataFrame)

    # quadratic FEM
    data_mfem2 = CSV.read("data/exp_pro/ref_rect_mfem2.csv",DataFrame)
    data_dfem2 = CSV.read("data/exp_pro/ref_rect_dfem2.csv",DataFrame)

    # linear XFEM
    data_mxfem1 = CSV.read("data/exp_pro/ref_rect_mxfem1.csv",DataFrame)
    data_dxfem1 = CSV.read("data/exp_pro/ref_rect_dxfem1.csv",DataFrame)

    # quadratic XFEM
    data_mxfem2 = CSV.read("data/exp_pro/ref_rect_mxfem2.csv",DataFrame)
    data_dxfem2 = CSV.read("data/exp_pro/ref_rect_dxfem2.csv",DataFrame)

    plot!(mass_plots, KRs, data1[!,iA], label="CutFEM p = $order",linecolor="#332288")#,markershape=:circle)
    plot!(mass_plots, KRs, data2[!,iA], label="AgFEM p = $order",linecolor="#CC6677")#,markershape=:cross)
    plot!(mass_plots, KRs, data3[!,iA], label="SBM p = $order",linecolor="#117733")#,markershape=:diamond)

    scatter!(mass_plots, data_me[!,1], data_me[!,2], label="experiments",markershape=:cross,markercolor=:black)#,markershape=:diamond)
    plot!(mass_plots,data_mfem1[!,1],data_mfem1[!,2],label="FEM p=1",linecolor="#44AA99", linestyle=:dash,markershape=:+,markercolor="#44AA99",markersize=msize)
    plot!(mass_plots,data_mfem2[!,1],data_mfem2[!,2],label="FEM p=2",linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)
    plot!(mass_plots,data_mxfem1[!,1],data_mxfem1[!,2],label="XFEM p=1",linecolor="#88CCEE", linestyle=:dash,markershape=:utriangle,markercolor="#88CCEE",markersize=msize)
    plot!(mass_plots,data_mxfem2[!,1],data_mxfem2[!,2],label="XFEM p=2",linecolor="#DDCC77", linestyle=:dash,markershape=:dtriangle,markercolor="#DDCC77",markersize=msize)
    plot!(mass_plots,data_mhpc[!,1],data_mhpc[!,2],label="HPC",linecolor="#882255", linestyle=:dash,markershape=:circle,markercolor="#882255",markersize=msize)
    

    scatter!(damping_plots, data_de[!,1], data_de[!,2], label="experiments",markershape=:cross,markercolor=:black)#,markershape=:diamond)
    plot!(damping_plots,data_dfem1[!,1],data_dfem1[!,2],label="FEM p=1",linecolor="#44AA99", linestyle=:dash,markershape=:+,markercolor="#44AA99",markersize=msize)
    plot!(damping_plots,data_dfem2[!,1],data_dfem2[!,2],label="FEM p=2",linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)
    plot!(damping_plots,data_dxfem1[!,1],data_dxfem1[!,2],label="XFEM p=1",linecolor="#88CCEE", linestyle=:dash,markershape=:utriangle,markercolor="#88CCEE",markersize=msize)
    plot!(damping_plots,data_dxfem2[!,1],data_dxfem2[!,2],label="XFEM p=2",linecolor="#DDCC77", linestyle=:dash,markershape=:dtriangle,markercolor="#DDCC77",markersize=msize)
    plot!(damping_plots,data_dhpc[!,1],data_dhpc[!,2],label="HPC",linecolor="#882255", linestyle=:dash,markershape=:circle,markercolor="#882255",markersize=msize)


    
    plot!(damping_plots, KRs, data1[!,iB], label="CutFEM p = $order",linecolor="#332288")#,markershape=:circle)
    plot!(damping_plots, KRs, data2[!,iB], label="AgFEM p = $order",linecolor="#CC6677")#,markershape=:cross)
    plot!(damping_plots, KRs, data3[!,iB], label="SBM p = $order",linecolor="#117733")#,markershape=:diamond)
    display(mass_plots)
    display(damping_plots)
    savefig(mass_plots,"plots/case41/A_$(name)_$order.png")
    savefig(damping_plots,"plots/case41/B_$(name)_$order.png")
end # function

# to = TimerOutput()
# # global variables      
# g = 9.81                    # [kg m/s²]: gravitational constant
# KRs = [0.1:0.025:2.0;]      # range of non-dimensional wave numbers
# order = 1                   # order of elements, either 1 or 2
# γg = 0.1                    # GP stabilization parameter 
# h = 0.0035                  # smallest element size in background mesh
# outputdir = "data/sims/"   # output directory

# # case specific variables
# R = 0.1                         # [m]: radius
# pmid = VectorValue(0.0,0.0)     # [m]: center point of radius
# ρV = π*R^2/2                    # [m]: area of a full horizontal cylinder (half domain)
# Ks = KRs./R                     # [m⁻¹]: range of wave numbers

# load background model, and geometry; then cut geometry into model
# model, geo = setup_domain(R,pmid,"data/meshes/background_shapes4.msh",Val(2),Val(:cylinder))
# model, geo = setup_domain(R,pmid,12.0,4.0,(1500,500),Val(2),Val(:cylinder))
# cutgeo, cutgeo_facets = cutting_model(model,geo)

# check d and n analytical
# degree = 2*order
# (Ω,Γ,_,_,_),(nΓ,),(dΩ,dΓ,dΓf,_,dΓo) = setup_interiors(model,cutgeo,degree)
# writevtk(Γ,"Gammasbm",cellfields=["d"=>d(pmid,R),"n"=>n(pmid)])

# run case for agfem, cutfem or sbm
# (aₐ,bₐ) = run_agfem(Ks, ρV, order, model, cutgeo, cutgeo_facets)
# (aₑ,bₑ) = run_cutfem(Ks, ρV, order, model, cutgeo, cutgeo_facets, γg, h)
# (aₛ,bₛ) = run_sbm(Ks, ρV, order, model, cutgeo, n(pmid), d(pmid,R))

# write_csv(aₐ,bₐ,outputdir*"agfem/cylHR0000_$order.csv";namex="A",namey="B")
# write_csv(aₑ,bₑ,outputdir*"cutfem/cylHR0000_$order.csv";namex="A",namey="B")
# write_csv(aₛ,bₛ,outputdir*"sbm/cylHR0000_$order.csv";namex="A",namey="B")

# plotter0000()
# show(to)



end # module