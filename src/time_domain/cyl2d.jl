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
using GridapGmsh
using LineSearches: BackTracking


# function definitions
function setup_domain(R::Float64,pmid::VectorValue)
    geo = disk(R, x0=pmid)
    model = GmshDiscreteModel("data/meshes/background_cyl.msh")
    model,geo
end # function

function cutting_model(model::DiscreteModel,geo)
    cutgeo = cut(model, !geo)
    cutgeo_facets = cut_facets(model, !geo)
    cutgeo, cutgeo_facets
end # function

function setup_spaces(order::Int64, model::DiscreteModel, Ω::Triangulation)
    reffeᵩ = ReferenceFE(lagrangian, Float64, order)
    reffeᵦ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    W = FESpace(Ω, reffeᵩ, vector_type=Vector{Float64})
    Φ = TransientTrialFESpace(W)
    V = ConstantFESpace(model; vector_type=Vector{Float64}, field_type=VectorValue{2, Float64})
    U = TransientTrialFESpace(V)
    E = FESpace(Ω, reffeᵦ, vector_type=Vector{Float64},dirichlet_tags=["seabed"])
    Η = TransientTrialFESpace(E,ηᵢ)
    W,Φ,V,U,E,Η
end # function

function setup_interiors(model,cutgeo,cutgeo_facets,degree)
    Ω = Interior(cutgeo, PHYSICAL)
    # writevtk(Ω,"model1")
    Ω⁻act = Interior(cutgeo, ACTIVE)
    Γ = EmbeddedBoundary(cutgeo)
    nΓ = get_normal_vector(Γ)
    Γf⁻ = BoundaryTriangulation(cutgeo_facets, tags=["surface"])
    dΓf⁻ = Measure(Γf⁻, degree)
    nΓf = get_normal_vector(Γf⁻)
    dΩ = Measure(Ω, degree)
    dΓ = Measure(Γ, degree)
    Γi = BoundaryTriangulation(model, tags=["inlet"])
    dΓi = Measure(Γi, degree)
    Γo = BoundaryTriangulation(model, tags=["outlet"])
    dΓo = Measure(Γo, degree)
    E⁰ = GhostSkeleton(cutgeo)
    dE⁰ = Measure(E⁰,degree)
    nE⁰ = get_normal_vector(E⁰)
    (Ω,Ω⁻act),(Γ,Γf⁻,Γi,Γo,E⁰),(nΓ,nΓf,nE⁰),(dΩ,dΓ,dΓf⁻,dΓi,dΓo,dE⁰)
end # function

R = 0.1524/2    # [m] radius; see Ito 1977 for experiment specifics
g = 9.81
order = 1
degree = 2*order


# initial conditions (zero displacement & velocity potential)
uᵢ(x,t) = VectorValue(0.0,R/3)               # [m] heave displacement Ito 1977
uᵢ(t::Real) = x -> uᵢ(x,t)
ϕᵢ(x,t) = 0.0               # [m²/s] velocity potential
ϕᵢ(t::Real) = x -> ϕᵢ(x,t)
ηᵢ(x,t) = VectorValue(0.0,0.0)               # [m] wave height
ηᵢ(t) = x -> ηᵢ(x,t)



model, geo = setup_domain(R,VectorValue(0.0,0.0))
cutgeo, cutgeo_facets = cutting_model(model,geo)
(Ω,Ω⁻act),(Γ,Γf⁻,Γi,Γo,E⁰),(n_Γ,n_Γf,nE⁰),(dΩ⁻,dΓ,dΓf⁻,dΓi,dΓo,dE⁰) = setup_interiors(model,cutgeo,cutgeo_facets,degree)
W,Φ,V,U,E,Η = setup_spaces(order, model, Ω⁻act)
X = TransientMultiFieldFESpace([Φ, U, Η])
Y = MultiFieldFESpace([W, V, E])

m_cyl = 0.03810023356932984#π*R^2 * ρᵦ/ρ /sum(∫(1.0)dΓ)    # [kg/m] mass

nz = VectorValue(0.0,1.0)

# Map from reference to deformed
# φ: Ω₀ ⟶ Ωₓ
id(x) = x
φ(t,u) = id + u

# Analytical solution evaluated from the reference configuration
# ϕ₀ = ϕₓ∘φ: Ω₀ ⟶ Ωₓ ⟶ ℛ
# ϕ₀(t,u) = ϕₓ(t)∘(φ(t,u))
# ϕ₀(t) = ϕ₀(t,u)

Ft(t,u) = ∇(φ(t,u))
invFt(t,u) = pinvJt∘(Ft(t,u))
nₓ(t,n) = u -> Gridap.Geometry.push_normal∘(invFt(t,u),n)
nₓ_Γ(t,u) = nₓ(t,n_Γ)(u)
nₓ_Γf(t,u) = nₓ(t,n_Γf)(u)
J(t,u) = meas∘(Ft(t,u))
C(t,u) = (ft->ft⋅ft')∘Ft(t,u)
J_Γ(t,u) = ((j,c,n)->j*sqrt(n⋅inv(c)⋅n))∘(J(t,u),C(t,u),n_Γ)
J_Γf(t,u) = ((j,c,n)->j*sqrt(n⋅inv(c)⋅n))∘(J(t,u),C(t,u),n_Γf)
∇ₓ(t,f,u) = invFt(t,u)⋅∇(f)
divₓ(t,f,u) = tr(∇ₓ(t,f,u))
∫ₓ(t,f) = u -> ∫(f*J(t,u))
∫ₓ_Γ(t,f) = u -> ∫(f*J_Γ(t,u))
∫ₓ_Γf(t,f) = u -> ∫(f*J_Γf(t,u))

a(t, (ϕ, u, η), (w, v, e)) = ∫(v ⋅ (∂tt(u) * m_cyl ) )dΓ 

b(t, (ϕ, u, η), (w, v, e)) = ∫( v ⋅ (∂t(ϕ) * n_Γ) )dΓ - 
                                ∫(w * ∂t(u) ⋅ n_Γ )dΓ + 
                                ∫(e ⋅ (∂t(ϕ) * n_Γf) )dΓf⁻ - 
                                ∫(w * ∂t(η) ⋅ n_Γf )dΓf⁻

c(t, (ϕ, u, η), (w, v, e)) = ∫(v ⋅ (( (g*nz) ⋅ u )*n_Γ))dΓ + 
                                ∫(∇(η)⊙∇(e) )dΩ⁻ +#-
                                # ∫(e⋅(∇(η)⋅n_Γ))dΓ - ∫(e⋅(∇(η)⋅n_Γf))dΓf⁻+ 
                                ∫(∇(ϕ)⋅∇(w) )dΩ⁻ + 
                                ∫((v ⋅ n_Γ)*0.5*(∇(ϕ)⋅∇(ϕ)))dΓ + 
                                ∫((e ⋅ n_Γf)*0.5*(∇(ϕ)⋅∇(ϕ)))dΓf⁻ + 
                                ∫(e ⋅ (( (g*nz) ⋅ η )*n_Γf))dΓf⁻

d(t, (ϕ, u, η), (w, v, e)) = ∫((w)*(-1.0)*(∇(ϕ)⋅∂t(u)))dΓ - 
                                ∫((w)*(∇(ϕ)⋅∂t(η)))dΓf⁻

rhs(t, (ϕ, u, η), (w, v, e)) = ∫(w * 0)dΩ⁻ + 
                                ∫(e ⋅ VectorValue(0.0,0.0))dΩ⁻

res(t, (ϕ, u, η), (w, v, e)) = a(t, (ϕ, u, η), (w, v, e)) + b(t, (ϕ, u, η), (w, v, e)) + c(t, (ϕ, u, η), (w, v, e)) + d(t, (ϕ, u, η), (w, v, e)) - rhs(t, (ϕ, u, η), (w, v, e))

jac(t, (ϕ, u, η), (dϕ, du, dη), (w, v, e)) = ∫(∇(dϕ)⋅∇(w) )dΩ⁻ + 
                                                ∫(∇(dη)⊙∇(e) )dΩ⁻+#-
                                                # ∫(e⋅(∇(dη)⋅n_Γ))dΓ +#- ∫(e⋅(∇(dη)⋅n_Γf))dΓf⁻+ 
                                                ∫( e ⋅ ( ((g*nz) ⋅ dη )*n_Γf))dΓf⁻ + 
                                                ∫( v ⋅ (( (g*nz) ⋅ du )*n_Γ))dΓ + 
                                                ∫((v ⋅ n_Γ)*0.5*(∇(dϕ)⋅∇(ϕ)))dΓ + 
                                                ∫( (e ⋅ n_Γf)*0.5*(∇(dϕ)⋅∇(ϕ)))dΓf⁻

jac_t(t, (ϕ, u, η), (dtϕ, dtu, dtη), (w, v, e)) = ∫(v ⋅ (dtϕ * n_Γ) )dΓ - 
                                                    ∫(w * dtu ⋅ n_Γ )dΓ + 
                                                    ∫(e ⋅ (dtϕ * n_Γf) )dΓf⁻ - 
                                                    ∫(w * dtη ⋅ n_Γf )dΓf⁻ +
                                                    ∫((w)*(-1.0)*(∇(ϕ)⋅dtu))dΓ - 
                                                    ∫((w)*(∇(ϕ)⋅dtη))dΓf⁻

jac_tt(t, (ϕ, u, η), (dttϕ, dttu, dttη), (w, v, e)) = ∫( v ⋅ (dttu * m_cyl ) )dΓ 


# mass_ql(t, (ϕ, u, η), (dttϕ, dttu, dttη), (w, v, e)) = ∫ₓ_Γ(t, v ⋅ (dttu * m_cyl ) )(η)dΓ 
# res_ql(t, (ϕ, u, η), (w, v, e)) = b(t, (ϕ, u, η), (w, v, e)) + c(t, (ϕ, u, η), (w, v, e)) - rhs(t, (w, v, e))

# time integrator variables
Δt = 1/100#10e-6#1/200                          # [s]: time step
t₀ = 0.0                            # [s]: zero time 
Tf = 5*Δt # 1.0#0.66 #3.0                           # [s]: final time
t_vec = [t₀:Δt:Tf]
ρ∞ = 0.0                    # time integrator constant corresponding to Newmark settings

# solver definition
ls = LUSolver()
nls = NLSolver(ls,show_trace=true, method=:newton, linesearch=BackTracking())

# op = TransientLinearFEOperator((c,b,a),rhs, X, Y;constant_forms=(true,true,true))
# op = TransientQuasilinearFEOperator(mass_ql, res_ql,(jac,jac_t,jac_tt),X,Y)
op = TransientFEOperator(res, (jac, jac_t, jac_tt), X,Y)

ode_solver = GeneralizedAlpha2(nls, Δt, ρ∞)

# interpolate initial conditions
x₀ = interpolate_everywhere([ϕᵢ(0.0), uᵢ(0.0), ηᵢ(0.0)], X(0.0))
v₀ = interpolate_everywhere([ϕᵢ(0.0), ηᵢ(0.0), ηᵢ(0.0)], X(0.0))

# solve
ϕhₜ = Gridap.solve(ode_solver, op, t₀, Tf, (x₀, v₀,v₀))

plot_u = Float64[]
push!(plot_u,R/3)
for (t,(ϕh,uh,ηh)) in ϕhₜ
    push!(plot_u,uh.free_values[1])
end


plt = plot(legend=:topright)
plot!(t_vec,plot_u./(R/3),xlabel="t[s]",ylabel="u[-]",title="Cylinder Heave Free Decay CutFEM")
display(plt)

stop

to = TimerOutput()
@timeit to "variables" begin
    # variables
    g = 9.81        # [kg/s²] gravitational constant
    L₁ = 6.0#18.0        # [m] length of domain
    d = 1.22        # [m] depth
    n₁ = 800#360        # [-] number of elements horizontally
    n₂ = 140         # [-] number of elements vertically
    R = 0.1524/2    # [m] radius; see Ito 1977 for experiment specifics
    println("Number of elements: ", n₁*n₂)
    ρ = 1000    # [kg/m³] density water
    ρᵦ = 500    # [kg/m³] density cylinder (half that of water)

    # Ghost penalty parameter
    γg = 0.1
    h = d/n₂

    # Damping
    μ₀ = 2.5
    x0 =  L₁/2*(0.75)
    Ld =  L₁/2-L₁/2*(0.75)
    μ₁(x::VectorValue) = μ₀*(1.0 - sin(π/2*(√(x[1]^2)-x0)/Ld))*(√(x[1]^2)<=x0)
    
    # initial conditions (zero displacement & velocity potential)
    uᵢ(x,t) = R/3               # [m] heave displacement Ito 1977
    uᵢ(t::Real) = x -> uᵢ(x,t)
    ϕᵢ(x,t) = 0.0               # [m²/s] velocity potential
    ϕᵢ(t::Real) = x -> ϕᵢ(x,t)
    ηᵢ(x,t) = 0.0               # [m] wave height
    ηᵢ(t) = x -> ηᵢ(x,t)
end



@timeit to "embedding" begin
    cutgeo = cut(model, (geo))
    cutgeo_facets = cut_facets(model, (geo))
    Ω = Interior(model) 
    Ω⁻ = Interior(cutgeo, PHYSICAL_OUT)
    Ω⁻act = Interior(cutgeo, ACTIVE_OUT)
    Γ = EmbeddedBoundary(cutgeo)
    n_Γ = -get_normal_vector(Γ)
    Γg = SkeletonTriangulation(cutgeo_facets, ACTIVE_OUT)
    n_Γg = get_normal_vector(Γg)
    Γg_i = GhostSkeleton(cutgeo, ACTIVE_OUT)
    n_Γg_i = get_normal_vector(Γg_i)
    
    Λ = GhostSkeleton(cutgeo, ACTIVE_OUT)
  
    # create integration space (triangulation) & Gauss quadratures (measure)
    order = 2
    degree = 2*order
    Γf⁻ = BoundaryTriangulation(cutgeo_facets, PHYSICAL_OUT, tags=["surface"])
    Γf_act = BoundaryTriangulation(cutgeo_facets, ACTIVE_OUT, tags=["surface"])
    Γf = BoundaryTriangulation(model, tags=["surface"])
    Γd = BoundaryTriangulation(model, tags=["damping"])
    dΓf_act = Measure(Γf_act, degree)
    dΓf⁻ = Measure(Γf⁻, degree)
    dΓd = Measure(Γd, degree)
    n_Γf = get_normal_vector(Γf⁻)
    dΩ⁻ = Measure(Ω⁻, degree)
    dΩ = Measure(Ω, degree)
    dΓ = Measure(Γ, degree)
    dΛ = Measure(Λ, degree)
    nΛ = get_normal_vector(Λ)

    # copied the degree of 8 from transientstokes - set to 2
    dΓg = Measure(Γg, degree)
    dΓg_i = Measure(Γg_i, degree)

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
    Wstd = FESpace(Ω⁻act, reffeᵩ)
    Φstd = TransientTrialFESpace(Wstd)
    Dstd = ConstantFESpace(model)
    Rstd = TransientTrialFESpace(Dstd)

    # final FE spaces
    X = TransientMultiFieldFESpace([Φstd, Rstd])
    Y = MultiFieldFESpace([Wstd, Dstd])
end

m_cyl = 0.03810023356932984#π*R^2 * ρᵦ/ρ /sum(∫(1.0)dΓ)    # [kg/m] mass

println("Initial displacement: ", uᵢ(VectorValue(0.0,0.0),0), " m")
println("Mass: ", m_cyl, " kg/m")

nz = VectorValue(0.0,1.0)

a(t, (ϕ, u), (w, v)) = ∫( v * (u * m_cyl ) )dΓ +
                        ∫(w * (ϕ/g ))dΓf⁻ +
                        ∫(w * μ₁ * (ϕ/g ))dΓd 

b(t, (ϕ, u), (w, v)) = ∫( v * (ϕ * (nz⋅n_Γ)) )dΓ - 
                        ∫( w * u * (nz⋅n_Γ) )dΓ 

c(t, (ϕ, u), (w, v)) = ∫( v * ( g * u * (nz⋅n_Γ) ))dΓ + 
                        ∫( ∇(ϕ)⋅∇(w) )dΩ⁻ +
                        ∫(γg*h^(2*2+1)*jump((∇∇(ϕ)⋅nΛ))⋅jump((∇∇(w)⋅nΛ)))dΛ +
                        ∫(γg*h^(2*1+1)*jump(∇(ϕ)⋅nΛ)*jump(∇(w)⋅nΛ))dΛ 

rhs(t, (w, v)) = ∫(w * 0)dΩ⁻

# time integrator variables
Δt = 1/50#10e-6#1/200                          # [s]: time step
t₀ = 0.0                            # [s]: zero time 
Tf = 2.0#0.66 #3.0                           # [s]: final time
t_vec = [t₀:Δt:Tf]
ρ∞ = 0.0                    # time integrator constant corresponding to Newmark settings

# solver definition
ls = LUSolver()
op = TransientLinearFEOperator((c,b,a),rhs, X, Y;constant_forms=(true,true,true))
ode_solver = GeneralizedAlpha2(ls, Δt, ρ∞)

# interpolate initial conditions
x₀ = interpolate_everywhere([ϕᵢ(0.0), uᵢ(0.0)], X(0.0))
v₀ = interpolate_everywhere([ϕᵢ(0.0), ηᵢ(0.0)], X(0.0))

# solve
ϕhₜ = Gridap.solve(ode_solver, op, t₀, Tf, (x₀, v₀,v₀))

folder="cylinder"
name="CutFEM"

pvd = createpvd("data/sims/"*folder*"/"*name)
pvd2 = createpvd("data/sims/"*folder*"/"*name*"_u")

plot_u = Float64[]
push!(plot_u,R/3)
@timeit to "solving" begin
    pvd[0] = createvtk(Ω⁻, "data/sims/"*folder*"/"*name*"_0.vtu", cellfields=[ "phih"=>ϕᵢ(0.0) ])
    pvd2[0] = createvtk(Γ, "data/sims/"*folder*"/"*name*"_u_0.vtu", cellfields=["uh"=>uᵢ(0.0)])
    data = DataFrame(Time=[0.0], Phi=[x₀[1].free_values], Disp=x₀[2].free_values)
    for (t,(ϕh,uh)) in ϕhₜ
        @timeit to "timesteps" begin
            global data
            pvd[t] = createvtk(Ω⁻,"data/sims/"*folder*"/"*name*"_$t"*".vtu",cellfields=["phih"=>ϕh])
            pvd2[t] = createvtk(Γ,"data/sims/"*folder*"/"*name*"_u_$t"*".vtu",cellfields=["uh"=>uh])
            data₁ = DataFrame(Time=[t], Phi=[interpolate_everywhere(ϕh,Wstd).free_values], Disp=[uh.free_values[1]])
            data = vcat(data,data₁)
            push!(plot_u,uh.free_values[1])
        end
    end
end
wsave("data/sims/"*folder*"/sol_cutfem_"*"$n₁"*"_"*"$n₂"*".jld2",Dict("df" => data))

plt = plot(legend=:topright)
plot!(t_vec,plot_u./(R/3),xlabel="t[s]",ylabel="u[-]",title="Cylinder Heave Free Decay CutFEM")
display(plt)
vtk_save(pvd)
vtk_save(pvd2)
show(to)
end