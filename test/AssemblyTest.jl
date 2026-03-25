using Test
using Gridap
using SparseArrays
using LinearAlgebra
using EmbeddedCoefficients
using EmbeddedCoefficients: _assemble_matrix, WeakForm, AssembledMatrices

# Setup: Create model, FEspaces and weak forms

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

W     = FESpace(model, ReferenceFE(lagrangian,Float64,1), vector_type=Vector{ComplexF64})
Φ     = TrialFESpace(W)
V = ConstantFESpace(model; vector_type=Vector{ComplexF64},
                        field_type=VectorValue{3,ComplexF64})
U = TrialFESpace(V)
spaces = FESpaces(W, Φ, V, U)

Ω = Triangulation(model)
dΩ = Measure(Ω, 2)

mock_vu = (u, v) -> ∫(u ⋅ v)dΩ
mock_vϕ = (u, v) -> ∫(u * v⋅VectorValue(1,0,0))dΩ
mock_wu = (u, v) -> ∫(u⋅VectorValue(1,0,0) * v)dΩ
mock_wϕ = (u, v) -> ∫(u * v)dΩ

# ============================================================================
# TEST SUITE
# ============================================================================

@testset "Matrix Assembly Tests" begin
    
    @testset "_assemble_matrix with Function" begin
        result = _assemble_matrix(mock_wϕ, Φ, W)
        @test result !== nothing
        @test isa(result, SparseMatrixCSC{ComplexF64, Int64})
        @test size(result) == (9,9)

        result = _assemble_matrix(mock_vu, U, V)
        @test result !== nothing
        @test isa(result, SparseMatrixCSC{ComplexF64, Int64})
        @test size(result) == (3,3)
    end
    
    @testset "_assemble_matrix with Nothing" begin
        result = _assemble_matrix(nothing, U, V)
        @test result === nothing
    end
    
    @testset "assemble_matrices with all functions present" begin
        a_wϕ = WeakForm(mock_wϕ, mock_wϕ)
        a_wu = WeakForm(mock_wu, mock_wu)
        a_vϕ = WeakForm(mock_vϕ, mock_vϕ)
        
        assembled = assemble_matrices(a_wϕ, a_wu, a_vϕ, spaces)
        
        # Check structure
        @test isa(assembled, AssembledMatrices)
        @test isa(assembled.dep, SystemMatrices)
        @test isa(assembled.indep, SystemMatrices)
        
        # Check dependent matrices
        @test assembled.dep.A_wϕ !== nothing
        @test assembled.dep.A_wu !== nothing
        @test assembled.dep.A_vϕ !== nothing
        
        # Check independent matrices
        @test assembled.indep.A_wϕ !== nothing
        @test assembled.indep.A_wu !== nothing
        @test assembled.indep.A_vϕ !== nothing
        
        # Check dimensions
        @test size(assembled.dep.A_wϕ) == (9, 9)
        @test size(assembled.dep.A_wu) == (9, 3)
        @test size(assembled.dep.A_vϕ) == (3, 9)
    end
    
    @testset "assemble_matrices with mixed functions and Nothing" begin
        a_wϕ = WeakForm(mock_wϕ, nothing)
        a_wu = WeakForm(nothing, mock_wu)
        a_vϕ = WeakForm(mock_vϕ, mock_vϕ)
        
        assembled = assemble_matrices(a_wϕ, a_wu, a_vϕ, spaces)
        
        # Check dependent part
        @test assembled.dep.A_wϕ !== nothing
        @test assembled.dep.A_wu === nothing  # Nothing in dep
        @test assembled.dep.A_vϕ !== nothing
        
        # Check independent part
        @test assembled.indep.A_wϕ === nothing  # Nothing in indep
        @test assembled.indep.A_wu !== nothing
        @test assembled.indep.A_vϕ !== nothing
    end
    
    @testset "assemble_matrices with all Nothing" begin
        a_wϕ = WeakForm(nothing, nothing)
        a_wu = WeakForm(nothing, nothing)
        a_vϕ = WeakForm(nothing, nothing)
        
        assembled = assemble_matrices(a_wϕ, a_wu, a_vϕ, spaces)
        
        # All matrices should be nothing
        @test assembled.dep.A_wϕ === nothing
        @test assembled.dep.A_wu === nothing
        @test assembled.dep.A_vϕ === nothing
        @test assembled.indep.A_wϕ === nothing
        @test assembled.indep.A_wu === nothing
        @test assembled.indep.A_vϕ === nothing
    end
    
    @testset "SystemMatrices construction" begin
        mat1 = sparse([1, 2], [1, 2], [1.0+0im, 2.0+0im], 2, 2)
        mat2 = sparse([1], [1], [1.0+0im], 2, 2)
        
        sys = SystemMatrices(mat1, mat2, nothing)
        
        @test sys.A_wϕ === mat1
        @test sys.A_wu === mat2
        @test sys.A_vϕ === nothing
    end
    
    @testset "AssembledMatrices structure" begin
        mat = sparse([1], [1], [1.0+0im], 2, 2)
        dep = SystemMatrices(mat, mat, nothing)
        indep = SystemMatrices(nothing, mat, mat)
        
        assembled = AssembledMatrices(dep, indep)
        
        @test assembled.dep === dep
        @test assembled.indep === indep
        @test assembled.dep.A_wϕ !== nothing
        @test assembled.indep.A_wϕ === nothing
    end
    
    @testset "Matrix sparsity and type correctness" begin
        a_wϕ = WeakForm(mock_wϕ, mock_wϕ)
        a_wu = WeakForm(mock_wu, mock_wu)
        a_vϕ = WeakForm(mock_vϕ, mock_vϕ)
        
        assembled = assemble_matrices(a_wϕ, a_wu, a_vϕ, spaces)
        
        # Check that returned matrices are actually sparse
        if assembled.dep.A_wϕ !== nothing
            @test issparse(assembled.dep.A_wϕ)
            @test eltype(assembled.dep.A_wϕ) == ComplexF64
        end
        
        if assembled.dep.A_wu !== nothing
            @test issparse(assembled.dep.A_wu)
            @test eltype(assembled.dep.A_wu) == ComplexF64
        end
        
        if assembled.dep.A_vϕ !== nothing
            @test issparse(assembled.dep.A_vϕ)
            @test eltype(assembled.dep.A_vϕ) == ComplexF64
        end
    end

end