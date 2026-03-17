struct FESpaces
    W::FESpace        # test space (radiation potential)
    Φ::FESpace        # trial space (radiation potential)
    V::FESpace        # test space (displacement: 3DOF in 2D or 6DOF in 3D)
    U::FESpace        # trial space (displacement: 3DOF in 2D or 6DOF in 3D)
end

function _num_hydro_dofs(model::DiscreteModel)
    N = num_cell_dims(model)
    N == 2 ? 3 : 6
end

function _setup_constant_space(model::DiscreteModel)
    V = ConstantFESpace(model; vector_type=Vector{ComplexF64},
                                field_type=VectorValue{_num_hydro_dofs(model),ComplexF64})
    U = TrialFESpace(V)
    V, U
end

function setup_spaces(order::Int64, model::DiscreteModel, Ω::Triangulation,
                        cutgeo::EmbeddedDiscretization, ::AGFEM,
                        threshold::Float64=1.0)
    strategy   = AggregateCutCellsByThreshold(threshold)
    aggregates = aggregate(strategy, cutgeo, cutgeo.geo, OUT)

    reffe = ReferenceFE(lagrangian, Float64, order)
    Wstd  = FESpace(Ω, reffe, vector_type=Vector{ComplexF64})
    W     = AgFEMSpace(Wstd, aggregates)
    Φ     = TrialFESpace(W)
    V, U  = _setup_constant_space(model)
    FESpaces(W, Φ, V, U)
end

function setup_spaces(order::Int64, model::DiscreteModel, Ω::Triangulation,
                        cutgeo::EmbeddedDiscretization,
                        ::Union{CUTFEM,SBM})
    reffe = ReferenceFE(lagrangian, Float64, order)
    W     = FESpace(Ω, reffe, vector_type=Vector{ComplexF64})
    Φ     = TrialFESpace(W)
    V, U  = _setup_constant_space(model)
    FESpaces(W, Φ, V, U)
end