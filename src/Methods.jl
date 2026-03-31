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

degree(m::EmbeddingMethod) = 2 * m.order

# Convenience constructors
AGFEM(; order=1, threshold=1.0)     = AGFEM(order, threshold)
CUTFEM(; order=1, γg=0.1, h=1.0)    = CUTFEM(order, γg, h)
SBM(; order=1)                      = SBM(order)

# Labels for plotting
label(m::AGFEM)  = "AgFEM (p=$(m.order))"
label(m::CUTFEM) = "CutFEM (p=$(m.order))"
label(m::SBM)    = "SBM (p=$(m.order))"

# Function to calculate h for CutFEM
function smallest_cell_size(model::DiscreteModel)
    trian = Triangulation(model)
    cell_coords = get_cell_coordinates(trian)
    
    # Define a function that computes minimum edge length for a single cell
    function cell_min_edge(coords::Vector) 
        npoints = length(coords)
        min_edge = Inf
        
        for i in 1:npoints
            for j in i+1:npoints
                dist = norm(coords[j] - coords[i])
                min_edge = min(min_edge, dist)
            end
        end
        
        return min_edge
    end
    
    # Apply the function to all cells using lazy evaluation, then collect results
    cell_sizes = lazy_map(cell_min_edge, cell_coords)
    min_h = minimum(collect(cell_sizes))
    
    return min_h
end