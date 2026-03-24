using STLCutters
using GridapEmbedded.Interfaces: SubFacetTriangulation
abstract type EmbeddedGeometry{N} end

# 2D
struct Circle <: EmbeddedGeometry{2}
    center::VectorValue{2,Float64}
    radius::Float64
end

struct Rectangle{N} <: EmbeddedGeometry{N}
    center::VectorValue{N,Float64}
    width::Float64    # B
    depth::Float64    # D
end

struct Triangle <: EmbeddedGeometry{2}
    center::VectorValue{2,Float64}
    depth::Float64    # D
    angle::Float64    # α in degrees
end

# 3D
struct Sphere <: EmbeddedGeometry{3}
    center::VectorValue{3,Float64}
    radius::Float64
end

struct VerticalCylinder <: EmbeddedGeometry{3}
    center::VectorValue{3,Float64}
    radius::Float64
    draft::Float64    # D
end

struct OC3 <: EmbeddedGeometry{3}
    folder::String
    center::VectorValue{3,Float64}
end

struct OC4 <: EmbeddedGeometry{3}
    folder::String
    center::VectorValue{3,Float64}
end

struct STL <: EmbeddedGeometry{3}
    folder::String
    center::VectorValue{3,Float64}
end

# Each geometry knows how to produce its Gridap levelset object
function gridap_geo(g::Circle)
    disk(g.radius, x0=g.center)
end

function gridap_geo(g::Rectangle)
    pmin = g.center - VectorValue(g.width/2, g.depth)
    quadrilateral(x0=pmin, d1=VectorValue(g.width, 0.0), d2=VectorValue(0.0, 2*g.depth))
end

function gridap_geo(g::Triangle)
    φ = g.angle * π / 180
    quadrilateral(x0=g.center, d1=VectorValue(g.depth*tan(φ), g.depth), d2=VectorValue(0.0, g.depth))
end

function gridap_geo(g::VerticalCylinder)
    bottom = g.center - VectorValue(0.0, 0.0, -g.draft)
    cylinder(g.radius; x0=bottom, v=VectorValue(0, 0, 1))
end

function gridap_geo(g::OC3)
    STLGeometry(g.folder)
end 

function gridap_geo(g::OC4)
    STLGeometry(g.folder)
end 

function gridap_geo(g::STL)
    STLGeometry(g.folder)
end 

# ===================================================
# Distances — levelset function
# ===================================================
function distances(bgmodel::DiscreteModel,
                    Γ::Triangulation,
                    geo::Union{Circle,Sphere}, degree::Int)

    Dspace   = num_point_dims(Γ)
    pmid     = geo.center
    QΓ       = CellQuadrature(Γ, degree)
    qcp      = get_cell_points(QΓ)
    z        = zero(VectorValue{Dspace, Float64})
    d_vec_cs = CellState(z, QΓ)
    n_vec_cs = CellState(z, QΓ)

    for (icell, cell) in enumerate(qcp.cell_phys_point)
        for (ipoint, point) in enumerate(cell)
            δ    = pmid - point          # vector from point to center
            absδ = sqrt(δ ⋅ δ)
            δ̂    = δ / absδ
            Xd   = pmid - geo.radius * δ̂     # point on boundary closest to x
            d    = Xd - point
            absd = sqrt(d ⋅ d)
            n    = absd > 0 ? d / absd : δ̂

            d_vec_cs.values[icell][ipoint] = d
            n_vec_cs.values[icell][ipoint] = n
        end
    end

    return d_vec_cs, n_vec_cs
end # Function

function _d(x,geo::Rectangle{2})
    dx = maximum([-geo.width/2 - x[1], 0.0, x[1] - geo.width/2])
    dy = maximum([-geo.depth - x[2], 0.0, x[2] - geo.depth])
    VectorValue(-dx,dy)
end # function

function _n(x,geo::Rectangle{2})
    dx = maximum([-geo.width/2 - x[1], 0.0, x[1] - geo.width/2])
    dy = maximum([-geo.depth - x[2], 0.0, x[2] - geo.depth])
    dist = √(dx^2+dy^2)
    VectorValue(-dx/dist,dy/dist)
end # function

function distances(bgmodel::DiscreteModel,
                    Γ::Triangulation,
                    geo::Rectangle{2}, degree::Int)

    Dspace   = num_point_dims(Γ)
    pmid     = geo.center
    QΓ       = CellQuadrature(Γ, degree)
    qcp      = get_cell_points(QΓ)
    z        = zero(VectorValue{Dspace, Float64})
    d_vec_cs = CellState(z, QΓ)
    n_vec_cs = CellState(z, QΓ)

    for (icell, cell) in enumerate(qcp.cell_phys_point)
        for (ipoint, point) in enumerate(cell)
            d_vec_cs.values[icell][ipoint] = _d(point,geo)
            n_vec_cs.values[icell][ipoint] = _n(point,geo)
        end
    end
    return d_vec_cs, n_vec_cs
end # function

# ===================================================
# Distances — STL
# ===================================================
# TO DO: verify and test STL cases
function distances(bgmodel::DiscreteModel,
                    Γ::BoundaryTriangulation,
                    geo::STLGeometry, degree::Int)

    topo           = Gridap.Geometry.get_grid_topology(bgmodel)
    D              = num_dims(topo)
    Dspace         = num_point_dims(Γ)
    QΓ             = CellQuadrature(Γ, degree)
    qcp            = get_cell_points(QΓ)

    cell_to_facets  = STLCutters.compute_cell_to_facets(bgmodel, STLCutters.get_stl(geo))
    face_to_cells   = Gridap.Geometry.get_faces(topo, D-1, D)
    face_to_facets  = STLCutters.compose_index_map(face_to_cells, cell_to_facets)
    face_to_facets  = face_to_facets[_face_to_bgface(Γ)]

    z        = zero(VectorValue{Dspace, Float64})
    d_vec_cs = CellState(z, QΓ)
    n_vec_cs = CellState(z, QΓ)

    for (icell, cell) in enumerate(qcp.cell_phys_point)
        Xc = STLCutters.closest_point(cell, geo,
                                        repeat(face_to_facets[icell], outer=length(cell)))
        for (ipoint, point) in enumerate(cell)
            δ = Xc[ipoint] - point
            d_vec_cs.values[icell][ipoint] = δ
            n_vec_cs.values[icell][ipoint] = δ / √(δ ⋅ δ)
        end
    end

    return d_vec_cs, n_vec_cs
end # function

# ===================================================
# Distance result container
# ===================================================
struct DistanceData
    d::CellState    # distance vector
    n::CellState    # normal vector
end # struct

# ===================================================
# Wrappers for each unfitted method: SBM 
# ===================================================
function compute_distances(::SBM, bgmodel::DiscreteModel,
                            geo::EmbeddedGeometry, Γ::BoundaryTriangulation,
                            degree::Int, t::Real)
    d, n = distances(bgmodel, Γ, geo, degree)
    return DistanceData(d, n)
end # function

# function compute_distances(::SBM, bgmodel::DiscreteModel,
#                             geo::EmbeddedGeometry, Γ::Triangulation,
#                             degree::Int, t::Real)
#     d, n = distances(bgmodel, Γ, geo, degree)
#     return DistanceData(d, n)
# end # function

function compute_distances(::SBM, bgmodel::DiscreteModel,
                            geo::STLGeometry, Γ::BoundaryTriangulation,
                            degree::Int, t::Real)
    d, n = distances(bgmodel, Γ, geo, degree)
    return DistanceData(d, n)
end # function

# AGFEM/CUTFEM — no distances needed
function compute_distances(::Union{AGFEM,CUTFEM}, bgmodel::DiscreteModel,
                            geo::EmbeddedGeometry, Γ::SubFacetTriangulation,
                            degree::Int, t::Real)
    return nothing
end