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

"""
Solve the system A_wϕ \\ A_wu for each column and project via A_vϕ.
Returns the y matrix of size (n_dof × n_dof).
"""
function solve_system(matrices::SystemMatrices)
    A_wϕ, A_wu, A_vϕ = matrices.A_wϕ, matrices.A_wu, matrices.A_vϕ
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