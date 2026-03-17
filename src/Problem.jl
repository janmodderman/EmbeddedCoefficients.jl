struct SimulationParams{N}
    domain::BackgroundMesh{N}
    geometry::EmbeddedGeometry{N}
    method::EmbeddingMethod
end

# # In problem.jl
# function setup_triangulations(p::SimulationParams, model, cutgeo, cutgeo_facets, degree)
#   setup_triangulations(model, cutgeo, cutgeo_facets, degree,
#                        p.method, p.domain.lateral_bc)
# end
