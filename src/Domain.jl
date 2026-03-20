abstract type BackgroundMesh{N} end
abstract type LateralTags end

# Cartesian background mesh (analytically defined)
struct CartesianDomain{N} <: BackgroundMesh{N}
    L₁::Float64
    L₂::Float64   # 3D only
    depth::Float64
    partition::Tuple
    lateral_tag::LateralTags
end

struct SymmetryInlet <: LateralTags end   # inlet = symmetry, outlet = wall
struct WallWall     <: LateralTags end    # both inlet and outlet = wall

# Convenience constructors
CartesianDomain2D(L₁, depth, partition; lateral_tag=WallWall()) = CartesianDomain{2}(L₁, 0.0, depth, partition, lateral_tag)
CartesianDomain3D(L₁, L₂, depth, partition; lateral_tag=WallWall()) = CartesianDomain{3}(L₁, L₂, depth, partition, lateral_tag)

# External Gmsh mesh
struct GmshDomain{N} <: BackgroundMesh{N}
    meshfile::String
end

function _build_cartesian_2d(d::CartesianDomain{2}, ::WallWall)
    pmin = Point(-d.L₁, -d.depth)
    pmax = Point(d.L₁,  0.0)
    CartesianDiscreteModel(pmin, pmax, d.partition)
end

function _build_cartesian_2d(d::CartesianDomain{2}, ::SymmetryInlet)
    pmin = Point(0.0,  -d.depth)
    pmax = Point(d.L₁,  0.0)
    CartesianDiscreteModel(pmin, pmax, (d.partition[1]÷2, d.partition[2]))
end

# Build the Gridap model + apply boundary tags
function setup_model(d::CartesianDomain{2})
    model = _build_cartesian_2d(d, d.lateral_tag)
    _tag_2d!(model, d.lateral_tag)
    return model
end

function setup_model(d::CartesianDomain{3})
    pmin  = Point(-d.L₁/2, -d.L₂/2, -d.depth)
    pmax  = Point(d.L₁/2,  d.L₂/2,  0.0)
    model = CartesianDiscreteModel(pmin, pmax, d.partition)
    _tag_3d!(model)
    return model
end

function setup_model(d::GmshDomain)
    GmshDiscreteModel(d.meshfile)   # tags come from the Gmsh file itself
end

# Private tagging helpers
function _tag_2d!(model, ::WallWall)
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "seabed",  [1,2,5])
    add_tag_from_tags!(labels, "walls",    [3,4,7,8])  # both lateral sides → wall
    add_tag_from_tags!(labels, "surface", [6])
end

function _tag_2d!(model, ::SymmetryInlet)
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "seabed",    [1,2,5])
    add_tag_from_tags!(labels, "symmetry",  [3,7])    # inlet side → symmetry
    add_tag_from_tags!(labels, "walls",      [4,8])    # outlet side → wall
    add_tag_from_tags!(labels, "surface",   [6])
end

function _tag_3d!(model)
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "seabed",  [21])
    add_tag_from_tags!(labels, "surface", [22])
    add_tag_from_tags!(labels, "walls", [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,23,24,25,26])
end

# container for results from cutting the background model
struct DiscreteCut 
    geo::Union{AnalyticalGeometry,STLGeometry}
    cut::EmbeddedDiscretization
    facets::EmbeddedFacetDiscretization
    stl::Union{STLEmbeddedDiscretization,Nothing}
end

function cut_model(model::DiscreteModel, geo::AnalyticalGeometry)
    return DiscreteCut(geo, cut(model, geo), cut_facets(model, geo), nothing)
end

function cut_model(model::DiscreteModel, geo::STLGeometry)
    cutgeo = cut(model, geo)
    return DiscreteCut(geo, cutgeo.cut, cut_facets(cutgeo), cutgeo)
end