abstract type EmbeddingMethod end 

struct AGFEM <: EmbeddingMethod
    order::Int64
    threshold::Float64
end

struct CUTFEM <: EmbeddingMethod
    order::Int64
    γg::Float64
    h::Float64
end

struct SBM <: EmbeddingMethod
    order::Int64
end

# Convenience constructors
AGFEM(; order=1, threshold=1.0)     = AGFEM(order, threshold)
CUTFEM(; order=1, γg=0.1, h=1.0)    = CUTFEM(order, γg, h)
SBM(; order=1)                      = SBM(order)

# Labels for plotting
label(m::AGFEM)  = "AgFEM (p=$(m.order))"
label(m::CUTFEM) = "CutFEM (p=$(m.order))"
label(m::SBM)    = "SBM (p=$(m.order))"