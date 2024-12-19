 module helpers 
 using Gridap, Gridap.Arrays, Gridap.FESpaces, Gridap.Fields, Gridap.Geometry, Gridap.TensorValues
 using GridapEmbedded, GridapEmbedded.Interfaces, GridapEmbedded.LevelSetCutters
 using GridapGmsh
 using WriteVTK
 using TimerOutputs
 using DrWatson
 using Plots
 using LinearAlgebra
 using LinearSolve
 using CSV, DataFrames

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
function setup_domain(B::Float64,D::Float64,pmin::VectorValue,mesh::String,::Val{2},::Val{:rectangle})
    model = GmshDiscreteModel(mesh)
    geo = quadrilateral(;x0=pmin,d1=VectorValue(B/2,0.0),d2=VectorValue(0.0,D))
    model,geo
end # function

function setup_domain(B::Float64,D::Float64,pmid::VectorValue,L₁::Float64,d::Float64,partition::Tuple,::Val{2},::Val{:rectangle})
    pmin = Point(0.0, -d)
    pmax = Point(L₁, 0.0)
    model = CartesianDiscreteModel(pmin, pmax, partition)
    geo = quadrilateral(;x0=pmid-VectorValue(B/2,D),d1=VectorValue(B,0.0),d2=VectorValue(0.0,D))
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

#==============================================================================================================================#

# CUTTING MODEL
function cutting_model(model,geo)
    cut(model, !geo), cut_facets(model, !geo)
end # function

#==============================================================================================================================#

# FE SPACES
# cell aggregation
function setup_spaces(order::Int64, model::DiscreteModel, Ω::Triangulation, cutgeo::EmbeddedDiscretization, dim::Int64)
    threshold = 1.0
    strategy = AggregateCutCellsByThreshold(threshold)
    aggregates = aggregate(strategy, cutgeo, cutgeo.geo)
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
    Ω = Interior(cutgeo, PHYSICAL)
    Ω⁻act = Interior(cutgeo, ACTIVE)
    Γ = EmbeddedBoundary(cutgeo)
    nΓ = get_normal_vector(Γ)
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
    Ω = Interior(cutgeo, IN)
    Γ = Interface(Interior(cutgeo,ACTIVE_OUT),Ω).⁻
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
J(u) = meas∘(Ft(∇(u)))
# J(u) = (Ft(∇(u)))⊙(Ft(∇(u)))

# pmid1(pmid,x) = pmid - x
# pmid1(pmid) = x -> pmid1(pmid,x)


# TRUE DISTANCE & NORMAL FUNCTIONS (SBM)
# horizontal cylinder (2D)
n(x,pmid) = (pmid-x)/√((pmid-x)⋅(pmid-x))
# n(x,pmid) = (pmid.-x)/√((pmid.-x)⋅(pmid.-x))
# n(x,pmid) = VectorValue(1.0,0.0)#(pmid.-x)/√((pmid.-x)⋅(pmid.-x))
d(x,pmid,R) =  (√((pmid-x)⋅(pmid-x))-R)*n(x,pmid)
# d(x,pmid,R) =  (√((pmid.-x)⋅(pmid.-x))-R)*n(x,pmid)
# d(x,pmid,R) = pmid(x).-x#(√((pmid.-x)⋅(pmid.-x))-R)*n(x,pmid)
n(pmid) = x -> n(x,pmid)
d(pmid,R) = x -> d(x,pmid,R)

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
    # y = zeros(ComplexF64,sz[1],sz[2])
    y1 = zeros(ComplexF64,2,2)
    for i in 1:sz[2]
        # u = Gridap.solve(LUSolver(),A_wϕ,A_wu[:,i])
        # y1[i,:] = A_vϕ*u 

        prob = LinearProblem(A_wϕ,A_wu[:,i])
        sol = LinearSolve.solve(prob)
        y1[i,:] = A_vϕ*sol.u 

        # global y[i,:] = sol.u         # To avoid copy this array of size n 
        #reshape(sol.u,(1,sz[1]))
    end # for
    # A_vϕ*y
    y1
end # function

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
function run_agfem(Ks, ρV, order, model, cutgeo, cutgeo_facets)
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

function run_cutfem(Ks, ρV, order, model, cutgeo, cutgeo_facets, γg, h;GPflag=true)
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

function run_sbm(Ks, ρV, order, model, cutgeo, n, d)
    added_mass = []
    added_damping = []
    degree = 2*order

    (Ω,_,_,_,_),(nΓ,),(dΩ,dΓ,dΓf,_,dΓo) = setup_interiors(model,cutgeo,degree)
    W,Φ,V,U = setup_spaces(order, model, Ω, num_dims(model))

    # required now to obtain gradient of d, TODO: find more better and elegant solution
    Vd= FESpace(Ω,ReferenceFE(lagrangian,VectorValue{num_dims(model),Float64},order)) # current implementation requires higher order to correctly get the gradient of the distance function
    dcf = interpolate_everywhere(CellField(d,Ω),Vd)
    # @show ∑(∫(∇(dcf)⊙∇(dcf))dΓ)
    # @show ∫(∇(dcf)⊙∇(dcf))dΓ
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


# y1 = Any[0.7185196587636757, 0.637885783829143, 0.5782845238528137, 0.5321282920210016, 0.4953140561041668, 0.46520883841893373, 0.4403467295260227, 0.41852520218703143, 0.39942378295282177, 0.385236168260485, 0.3727114543606869, 0.3588885419250259, 0.3498156540606156, 0.34274609545264795, 0.3349563100718828, 0.323928085443556, 0.3172377054643731, 0.3176846147773164, 0.3174754047370269, 0.3090674990860375, 0.30825406843078573, 0.3106380147651668, 0.30031216904500363, 0.30178451920041666, 0.2995910335629596, 0.30239694485110036, 0.29415860427140017, 0.2977771660024027, 0.30063556559227045, 0.29322934392486977, 0.29723969507722137, 0.29259773813367335, 0.30473745767585464, 0.30182437694209424, 0.3056905653110921, 0.3073123567567107, 0.29917528718050246, 0.28599360154185405, 0.3049216287102037, 0.3136525865277211, 0.3141884024998921, 0.31060642837562286, 0.31052052221371806, 0.3023322430439674, 0.3093107838247952, 0.3172835351747758, 0.3204271141613397, 0.3252692443402839, 0.31880333938428623, 0.30372109820154564, 0.3270818556626764, 0.33496043652585483, 0.3170704885415966, 0.33823673749981753, 0.32164216374134647, 0.3234400248574641, 0.33756613740467617, 0.3433402920936833, 0.33254072611903696, 0.3303335774855977, 0.3332229725406265, 0.35457922582617457, 0.33967242894886, 0.34275909735307936, 0.32908888447649914, 0.365383884356667, 0.3562199762683824, 0.3302835407394619, 0.3523576209605079, 0.3715440755695462, 0.34636279154917865, 0.35175769637181326, 0.3515554109444105, 0.3763505304602129, 0.3500799935075541, 0.3546296148364313, 0.3708721220723904]
# y2 = Any[0.8591286608534882, 0.8100375816149905, 0.7646730887642254, 0.7236837379279506, 0.6868019313753245, 0.6528136807890572, 0.6216243444548182, 0.5921228046562694, 0.5663605594843102, 0.5428865916721924, 0.5173097207858156, 0.49630679855468474, 0.4772868741564141, 0.45712300123813915, 0.4366214859497943, 0.4189085351925913, 0.4092200375884153, 0.3949941067253721, 0.3765949350434006, 0.3607037976484927, 0.3515337670084326, 0.33449015952559535, 0.31823944891342715, 0.31367866579639153, 0.30012704569520515, 0.2875185064198466, 0.2774130896211255, 0.27128631343643156, 0.2607014266435948, 0.24693453124022188, 0.24512019958876555, 0.23838240180696427, 0.23130423426577704, 0.22031087294535545, 0.2123778812137151, 0.19823762768580638, 0.1872619613275965, 0.19109789764909385, 0.20182703180984046, 0.18182927550398517, 0.17382133173132125, 0.16148280420681074, 0.15998619203483633, 0.1517657133386099, 0.15776804676516304, 0.15850675901159111, 0.1444679881853747, 0.13391343287332114, 0.1241933694609342, 0.13287118454133787, 0.1398664328239154, 0.11651487706845237, 0.12039288234199447, 0.11444310724427467, 0.10171379097088533, 0.12061653314492866, 0.10946349289550752, 0.10134708618942319, 0.08988294265019924, 0.09632267646721837, 0.10902602481342509, 0.08502584815464123, 0.08685486112569402, 0.07850670485580602, 0.09036676159042172, 0.09531550091825175, 0.0616363704816172, 0.06954404108220237, 0.09701476454897143, 0.06854937498351144, 0.05832509916009539, 0.0698808469380804, 0.08553909942500625, 0.05630157320847624, 0.04922834201341768, 0.07925727793164851, 0.06590800274309089]

# y1 = Any[0.7182944737823116, 0.6373841879157774, 0.5778316434856834, 0.5316972682296798, 0.4947173272349353, 0.4644352226143809, 0.43913650828825596, 0.4178973078233463, 0.39985213185589175, 0.38447415303437216, 0.37115840047690357, 0.35967309659053026, 0.3495140493779759, 0.3409093561047104, 0.3337513666616771, 0.32724222337928327, 0.32204537742556855, 0.31697190823720856, 0.31233095544183204, 0.30917680265277403, 0.3067933825801121, 0.30439793268234266, 0.3021376919748066, 0.3006578737050532, 0.2987016973330239, 0.2976269498867228, 0.29729214610217847, 0.29639068647548683, 0.29804229723854453, 0.29947470525807046, 0.29830252349490227, 0.29669975870707227, 0.2967536285330529, 0.29909227656687576, 0.2984836330125522, 0.2998086343983736, 0.30017976892997567, 0.3044838104773701, 0.3062518517537721, 0.3052825857144318, 0.3063641801521448, 0.3085740373415959, 0.3090744896868355, 0.30723750779130066, 0.31013116250258504, 0.3146376260068354, 0.3179348678792444, 0.31792985431858806, 0.3168786429073354, 0.3203662303433941, 0.3227540967950564, 0.3216129551772921, 0.32399439053090456, 0.32418401941837194, 0.3288846125358081, 0.33470567146298885, 0.3295477292395301, 0.33121285148370494, 0.33771679312115016, 0.33398016818501935, 0.3337449709745275, 0.3435698505620953, 0.3397101399133165, 0.34638963475752976, 0.3437902059645227, 0.34390154898689346, 0.3468015493285692, 0.3447848888471356, 0.35697943418012, 0.35107597830502735, 0.35366846452441797, 0.35543148086289367, 0.3509938770732782, 0.3614196542878078, 0.35811704746833456, 0.36605013476629294, 0.3578027460532878]
# y2 = Any[0.8586628026265243, 0.8100351788804363, 0.7647385907826654, 0.723771283286817, 0.6867183558593003, 0.652910995703574, 0.621980728986709, 0.5934298573736749, 0.566928069721861, 0.5421965000506397, 0.5190143676043351, 0.49723372735178417, 0.47685374825698906, 0.4579419524244054, 0.43970402539882764, 0.4225362578123775, 0.40609090672617376, 0.3901151926530417, 0.37577862357560576, 0.36244303769763614, 0.34883057648663807, 0.33567216477263984, 0.32356949727210055, 0.31183581582747993, 0.30024726349083236, 0.29066288118408906, 0.28047175318657686, 0.2712265059093779, 0.26317727170539507, 0.25200168584581334, 0.24175989731697992, 0.2330100099516471, 0.228267565112033, 0.21951067113604683, 0.2120025830133877, 0.20669589961018162, 0.1999986297035394, 0.19536704834844118, 0.185306506020133, 0.17892314086729874, 0.17416516867446855, 0.16848421605065322, 0.16108262140297097, 0.1578260446310207, 0.15643452718352227, 0.1511704755008185, 0.14540591487751348, 0.1378587849826561, 0.13415914416710287, 0.13267836954831158, 0.12626151512083983, 0.12234727646072033, 0.11998653187322103, 0.118487721702337, 0.11625244603550335, 0.10892587260321825, 0.1023322367571235, 0.1067914968128559, 0.09904993414682163, 0.0947798524854675, 0.09787233702711859, 0.09375017686571116, 0.0898117311689404, 0.08693959818822779, 0.08124394290390233, 0.08336026237970534, 0.07779790062310567, 0.08264839245607365, 0.0772692772306618, 0.07036884672456134, 0.07236229592344588, 0.06585534662722622, 0.07041488173764704, 0.0684877607939143, 0.06582843955808401, 0.06107930176961239, 0.05633169978122957]

# 

# TESTING
function plotter0000()
KRs = [0.1:0.025:2.0;]      # range of non-dimensional wave numbers

mass_plots = plot(title="Added Mass vs KR", xlabel="KR [-]", ylabel="Added Mass Mₐ [-]")
damping_plots = plot(title="Added Damping vs KR", xlabel="KR [-]", ylabel="Added Damping Cₐ [-]")
data1 = CSV.read("data/sims/cutfem/cylHR0000_1.csv",DataFrame)
data2 = CSV.read("data/sims/agfem/cylHR0000_1.csv",DataFrame)
data3 = CSV.read("data/sims/sbm/cylHR0000_1.csv",DataFrame)
# 
iA = 4
iB = 8
msize = 4
plot!(mass_plots, KRs, data1[!,iA], label="CutFEM",linecolor="#332288")#,markershape=:circle)
plot!(mass_plots, KRs, data2[!,iA], label="AgFEM",linecolor="#CC6677")#,markershape=:cross)
plot!(mass_plots, KRs, data3[!,iA], label="SBM",linecolor="#117733")#,markershape=:diamond)
refdatam = CSV.read("data/exp_pro/ref_cyl_me_HR0000.csv",DataFrame)
refdatad = CSV.read("data/exp_pro/ref_cyl_de_HR0000.csv",DataFrame)

plot!(mass_plots,refdatam[!,1],refdatam[!,2],label="h/R=0.000",linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)

plot!(damping_plots, KRs, data1[!,iB], label="CutFEM",linecolor="#332288")#,markershape=:circle)
plot!(damping_plots, KRs, data2[!,iB], label="AgFEM",linecolor="#CC6677")#,markershape=:cross)
plot!(damping_plots, KRs, data3[!,iB], label="SBM",linecolor="#117733")#,markershape=:diamond)
plot!(damping_plots,refdatad[!,1],refdatad[!,2],label="h/R=0.000",linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)
display(mass_plots)
display(damping_plots)
end # function

to = TimerOutput()
# global variables      
g = 9.81                    # [kg m/s²]: gravitational constant
KRs = [0.1:0.025:2.0;]      # range of non-dimensional wave numbers
order = 1                   # order of elements, either 1 or 2
γg = 0.1                    # GP stabilization parameter 
h = 0.0035                  # smallest element size in background mesh
outputdir = "data/sims/"   # output directory

# case specific variables
R = 0.1                         # [m]: radius
pmid = VectorValue(0.0,0.0)     # [m]: center point of radius
ρV = π*R^2/2                    # [m]: area of a full horizontal cylinder (half domain)
Ks = KRs./R                     # [m⁻¹]: range of wave numbers

# load background model, and geometry; then cut geometry into model
# model, geo = setup_domain(R,pmid,"data/meshes/background_shapes4.msh",Val(2),Val(:cylinder))
model, geo = setup_domain(R,pmid,12.0,4.0,(1500,500),Val(2),Val(:cylinder))
cutgeo, cutgeo_facets = cutting_model(model,geo)

# check d and n analytical
# degree = 2*order
# (Ω,Γ,_,_,_),(nΓ,),(dΩ,dΓ,dΓf,_,dΓo) = setup_interiors(model,cutgeo,degree)
# writevtk(Γ,"Gammasbm",cellfields=["d"=>d(pmid,R),"n"=>n(pmid)])

# run case for agfem, cutfem or sbm
# (aₐ,bₐ) = run_agfem(Ks, ρV, order, model, cutgeo, cutgeo_facets)
# (aₑ,bₑ) = run_cutfem(Ks, ρV, order, model, cutgeo, cutgeo_facets, γg, h)
(aₛ,bₛ) = run_sbm(Ks, ρV, order, model, cutgeo, n(pmid), d(pmid,R))

# write_csv(aₐ,bₐ,outputdir*"agfem/cylHR0000_$order.csv";namex="A",namey="B")
# write_csv(aₑ,bₑ,outputdir*"cutfem/cylHR0000_$order.csv";namex="A",namey="B")
write_csv(aₛ,bₛ,outputdir*"sbm/cylHR0000_$order.csv";namex="A",namey="B")

plotter0000()
# show(to)



end # module