struct SystemMatrices
    A_wϕ::Union{SparseMatrixCSC{ComplexF64, Int64}, Nothing}
    A_wu::Union{SparseMatrixCSC{ComplexF64, Int64}, Nothing}
    A_vϕ::Union{SparseMatrixCSC{ComplexF64, Int64}, Nothing}
end

struct AssembledMatrices
    dep::SystemMatrices
    indep::SystemMatrices
end

function _assemble_matrix(fun::Function, U::FESpace , V::FESpace )
    return assemble_matrix(fun, U, V)
end

function _assemble_matrix(fun::Nothing, U::FESpace , V::FESpace )
    return nothing
end

function assemble_matrices(a_wϕ::WeakForm, a_wu::WeakForm, a_vϕ::WeakForm, spaces::FESpaces)
    return AssembledMatrices(SystemMatrices(_assemble_matrix(a_wϕ.dep, spaces.Φ, spaces.W),
                                            _assemble_matrix(a_wu.dep, spaces.U, spaces.W),
                                            _assemble_matrix(a_vϕ.dep, spaces.Φ, spaces.V)), 
                            SystemMatrices(_assemble_matrix(a_wϕ.indep, spaces.Φ, spaces.W),
                                            _assemble_matrix(a_wu.indep, spaces.U, spaces.W),
                                            _assemble_matrix(a_vϕ.indep, spaces.Φ, spaces.V)))
end