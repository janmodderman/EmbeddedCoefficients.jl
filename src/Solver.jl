struct SimulationParams{N}
    domain::BackgroundMesh{N}
    geometry::EmbeddedGeometry{N}
    method::EmbeddingMethod
end

struct SimulationState
    meas::Measures
    norms::Normals
    spaces::FESpaces
    dist::Union{DistanceData,Nothing}
end

struct SolverCache
    n::Int64        
    y::Matrix{<:ComplexF64}
    u::Vector{<:ComplexF64}
    ns::Gridap.Algebra.LUNumericalSetup
end

function setup_simulation(params::SimulationParams)
    # 1. Domain
    model               = setup_model(params.domain)
    cut                 = cut_model(model, gridap_geo(params.geometry))

    # 2. Triangulations
    tris, meas, norms   = setup_triangulations(model, cut.cut, cut.facets, 
                                                degree(params.method), params.method, 
                                                params.domain.lateral_tag)

    # 3. Spaces
    spaces              = setup_spaces(params.method.order, model, 
                                        tris.Ω_act, cut, params.method)

    # 4. Distance functions
    dist                = compute_distances(params.method, model, 
                                            params.geometry, tris.Γ, degree(params.method))
    return SimulationState(meas, norms, spaces, dist)
end     

function pre_assemble(params::SimulationParams, setup::SimulationState)
    # 5. Weak forms
    a_wϕ = make_a_wϕ(params.method, setup.meas, setup.norms, setup.dist)
    a_wu = make_a_wu(params.method, params.geometry, setup.norms, setup.meas, setup.dist)
    a_vϕ = make_a_vϕ(params.method, params.geometry, setup.norms, setup.meas, setup.dist)

    # 6. Assembly
    return assemble_matrices(a_wϕ, a_wu, a_vϕ, setup.spaces)
end 

function setup_cache(matrices::AssembledMatrices)
    A_wϕ    = matrices.dep.A_wϕ + matrices.indep.A_wϕ
    ss      = symbolic_setup(LUSolver(), A_wϕ)
    ns      = numerical_setup(ss, A_wϕ)
    A_wu    = matrices.dep.A_wu
    n       = size(A_wu)[2]
    u       = similar(Vector{ComplexF64}(A_wu[:,1]))
    y       = zeros(ComplexF64, n, n)
    return SolverCache(n,y,u,ns)
end

function solve_system(cache::SolverCache, matrices::AssembledMatrices, k::Float64, ω::Float64)
    A_wϕ, A_wu, A_vϕ = matrices.dep.A_wϕ*k + matrices.indep.A_wϕ, matrices.dep.A_wu*ω, matrices.dep.A_vϕ*ω
    numerical_setup!(cache.ns, A_wϕ)
    for i in 1:cache.n
        Gridap.Algebra.solve!(cache.u, cache.ns, Vector{ComplexF64}(A_wu[:,i]))
        cache.y[i,:] = A_vϕ * cache.u
    end
    return cache.y
end

# TODO: look into variable domain size per wavenumber k
# TODO: do we want to pass ρV or calculate it inside?
function coeff_solve(cache::SolverCache, matrices::AssembledMatrices, ω::Float64, k::Float64, ρV::Float64)
    y = solve_system(cache, matrices, k, ω)
    return hydro_coeffs(y, ω, ρV)
end