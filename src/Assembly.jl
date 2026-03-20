struct SystemMatrices
    A_wϕ
    A_wu
    A_vϕ
end

function assemble_matrices(a_wϕ::Function, a_wu::Function, a_vϕ::Function,
                            spaces::FESpaces)
    A_wϕ = assemble_matrix(a_wϕ, spaces.Φ, spaces.W)
    A_wu = assemble_matrix(a_wu, spaces.U, spaces.W)
    A_vϕ = assemble_matrix(a_vϕ, spaces.Φ, spaces.V)
    SystemMatrices(A_wϕ, A_wu, A_vϕ)
end

struct AssembledMatrices
    dep::SystemMatrices
    indep::SystemMatrices
end

function _assemble_matrix(fun::Function, U, V)
    assemble_matrix(fun,U,V)
end

function _assemble_matrix(fun::Nothing, U, V)
    spzeros(ComplexF64, num_free_dofs(V), num_free_dofs(U))
end

function assemble_matrices(a_wϕ::WeakForm, a_wu::WeakForm, a_vϕ::WeakForm, spaces::FESpaces)
    AssembledMatrices(SystemMatrices(_assemble_matrix(a_wϕ.dep, spaces.Φ, spaces.W),
                            _assemble_matrix(a_wu.dep, spaces.U, spaces.W),
                            _assemble_matrix(a_vϕ.dep, spaces.Φ, spaces.V)), 
            SystemMatrices(_assemble_matrix(a_wϕ.indep, spaces.Φ, spaces.W),
                            _assemble_matrix(a_wu.indep, spaces.U, spaces.W),
                            _assemble_matrix(a_vϕ.indep, spaces.Φ, spaces.V)))
end

"""
Solve the system A_wϕ \\ A_wu for each column and project via A_vϕ.
Returns the y matrix of size (n_dof × n_dof).
"""
function solve_system(matrices::AssembledMatrices, k::Float64, ω::Float64)

    # A_wϕ, A_wu, A_vϕ = matrices.dep.A_wϕ*k + matrices.indep.A_wϕ, matrices.dep.A_wu*k, matrices.dep.A_vϕ*k
    A_wϕ, A_wu, A_vϕ = matrices.dep.A_wϕ*k + matrices.indep.A_wϕ, matrices.dep.A_wu*ω + matrices.indep.A_wu, matrices.dep.A_vϕ*ω + matrices.indep.A_vϕ

    sz = size(A_wu)
    y  = zeros(ComplexF64, sz[2], sz[2])
    ss = symbolic_setup(LUSolver(), A_wϕ)
    ns = numerical_setup(ss, A_wϕ)
    u  = similar(Vector{ComplexF64}(A_wu[:,1]))
    for i in 1:sz[2]
        u       = Gridap.Algebra.solve!(u, ns, Vector{ComplexF64}(A_wu[:,i]))
        y[i,:]  = A_vϕ * u
    end
    y
end