using STLCutters
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

# ===================================================
# Distances — levelset function
# ===================================================
function distances(bgmodel::DiscreteModel,
                    Γ::BoundaryTriangulation,
                    geo::Union{Circle,Sphere}, degree::Int)

    Dspace   = num_point_dims(Γ)
    pmid     = geo.center
    QΓ       = CellQuadrature(Γ, degree)
    qcp      = get_cell_points(QΓ)
    z        = zero(VectorValue{Dspace, Float64})
    d_vec_cs = CellState(z, QΓ)
    n_vec_cs = CellState(z, QΓ)
    Xd_cs    = CellState(z, QΓ)

    for (icell, cell) in enumerate(qcp.cell_phys_point)
        for (ipoint, point) in enumerate(cell)
            δ    = pmid - point          # vector from point to center
            absδ = sqrt(δ ⋅ δ)
            δ̂    = δ / absδ
            Xd   = pmid - geo.R * δ̂     # point on boundary closest to x
            d    = Xd - point
            absd = sqrt(d ⋅ d)
            n    = absd > 0 ? d / absd : δ̂

            d_vec_cs.values[icell][ipoint] = d
            n_vec_cs.values[icell][ipoint] = n
            Xd_cs.values[icell][ipoint]    = Xd
        end
    end

    return d_vec_cs, n_vec_cs, Xd_cs
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
    Xd_cs = CellState(z, QΓ)

    for (icell, cell) in enumerate(qcp.cell_phys_point)
        Xc = STLCutters.closest_point(cell, geo,
                                        repeat(face_to_facets[icell], outer=length(cell)))
        for (ipoint, point) in enumerate(cell)
            δ = Xc[ipoint] - point
            d_vec_cs.values[icell][ipoint] = δ
            n_vec_cs.values[icell][ipoint] = δ / √(δ ⋅ δ)
            Xd_cs.values[icell][ipoint] = point + d_vec_cs.values[icell][ipoint]
        end
    end

    return d_vec_cs, n_vec_cs, Xd_cs
end # function

# ===================================================
# Distance result container
# ===================================================
struct DistanceData
    d::CellState    # distance vector
    n::CellState    # normal vector
    Xd::CellState   # shifted point
end # struct

# ===================================================
# Wrappers for each unfitted method: SBM or WSBM
# ===================================================
function compute_distances(::SBM, bgmodel::DiscreteModel,
                            geo::EmbeddedGeometry, Γ::BoundaryTriangulation,
                            fun::Function, degree::Int, t::Real)
    d, n, Xd = distances(bgmodel, Γ, geo, degree)
    return DistanceData(d, n, Xd)
end # function

function compute_distances(::SBM, bgmodel::DiscreteModel,
                            geo::STLGeometry, Γ::BoundaryTriangulation,
                            fun::Function, degree::Int, t::Real)
    d, n, Xd = distances(bgmodel, Γ, geo, degree)
    return DistanceData(d, n, Xd)
end # function