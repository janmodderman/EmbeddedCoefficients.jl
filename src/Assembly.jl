struct SystemMatrices
    A_wϕ
    A_wu
    A_vϕ
end

struct AssembledMatrices
    dep::SystemMatrices
    indep::SystemMatrices
end

function _assemble_matrix(fun::Function, U, V)
    assemble_matrix(fun,U,V)
end

function _assemble_matrix(fun::Nothing, U, V)
    nothing
end

function assemble_matrices(a_wϕ::WeakForm, a_wu::WeakForm, a_vϕ::WeakForm, spaces::FESpaces)
    AssembledMatrices(SystemMatrices(_assemble_matrix(a_wϕ.dep, spaces.Φ, spaces.W),
                            _assemble_matrix(a_wu.dep, spaces.U, spaces.W),
                            _assemble_matrix(a_vϕ.dep, spaces.Φ, spaces.V)), 
            SystemMatrices(_assemble_matrix(a_wϕ.indep, spaces.Φ, spaces.W),
                            _assemble_matrix(a_wu.indep, spaces.U, spaces.W),
                            _assemble_matrix(a_vϕ.indep, spaces.Φ, spaces.V)))
end

struct SolverCache
    n::Int64        
    y::Matrix{<:ComplexF64}
    u::Vector{<:ComplexF64}
    ns::Gridap.Algebra.LUNumericalSetup
end

function setup_cache(matrices::AssembledMatrices)
    A_wϕ = matrices.dep.A_wϕ + matrices.indep.A_wϕ
    ss = symbolic_setup(LUSolver(), A_wϕ)
    ns = numerical_setup(ss, A_wϕ)
    A_wu = matrices.dep.A_wu
    n = size(A_wu)[2]
    u  = similar(Vector{ComplexF64}(A_wu[:,1]))
    y  = zeros(ComplexF64, n, n)
    SolverCache(n,y,u,ns)
end

function solve_system(cache::SolverCache, matrices::AssembledMatrices, k::Float64, ω::Float64)
    A_wϕ, A_wu, A_vϕ = matrices.dep.A_wϕ*k + matrices.indep.A_wϕ, matrices.dep.A_wu*ω, matrices.dep.A_vϕ*ω
    numerical_setup!(cache.ns, A_wϕ)
    for i in 1:cache.n
        Gridap.Algebra.solve!(cache.u, cache.ns, Vector{ComplexF64}(A_wu[:,i]))
        cache.y[i,:]  = A_vϕ * cache.u
    end
    cache.y
end