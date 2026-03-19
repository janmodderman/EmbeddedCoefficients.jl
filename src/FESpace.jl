struct FESpaces
    W::FESpace        # test space (radiation potential)
    Φ::FESpace        # trial space (radiation potential)
    V::FESpace        # test space (displacement: 3DOF in 2D or 6DOF in 3D)
    U::FESpace        # trial space (displacement: 3DOF in 2D or 6DOF in 3D)
end

# Returns degrees of freedom of the structure for dimension size of the problem: 2D => 3DOF, 3D => 6DOF
function _num_hydro_dofs(model::DiscreteModel)
    N = num_cell_dims(model)
    return N == 2 ? 3 : 6
end

# Returns test and trial FE Space with a complex valued single unknown of size _num_hydro_dofs(model)
function _setup_constant_space(model::DiscreteModel)
    V = ConstantFESpace(model; vector_type=Vector{ComplexF64},
                                field_type=VectorValue{_num_hydro_dofs(model),ComplexF64})
    U = TrialFESpace(V)
    return V, U
end

# Returns Aggregated test and trial FE Space for the velocity potential, and test and trial space for displacement
function setup_spaces(order::Int64, model::DiscreteModel, Ω::Triangulation,
                        cutgeo::Union{EmbeddedDiscretization}, ::AGFEM,
                        threshold::Float64=1.0)
    strategy   = AggregateCutCellsByThreshold(threshold)
    aggregates = aggregate(strategy, cutgeo, cutgeo.geo, OUT)

    reffe = ReferenceFE(lagrangian, Float64, order)
    Wstd  = FESpace(Ω, reffe, vector_type=Vector{ComplexF64})
    W     = AgFEMSpace(Wstd, aggregates)
    Φ     = TrialFESpace(W)
    V, U  = _setup_constant_space(model)
    return FESpaces(W, Φ, V, U)
end

# Returns test and trial FE Space for the velocity potential, and test and trial space for displacement
function setup_spaces(order::Int64, model::DiscreteModel, Ω::Triangulation,
                        cutgeo::EmbeddedDiscretization,
                        ::Union{CUTFEM,SBM})
    reffe = ReferenceFE(lagrangian, Float64, order)
    W     = FESpace(Ω, reffe, vector_type=Vector{ComplexF64})
    Φ     = TrialFESpace(W)
    V, U  = _setup_constant_space(model)
    return FESpaces(W, Φ, V, U)
end