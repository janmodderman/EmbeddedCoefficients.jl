# TODO: use alternative for type NamedTuple? Difficulty speficying for multiple dispatch names.
struct Triangulations
    Ω::Triangulation
    Ω_act::Triangulation
    Γ::Union{SubFacetTriangulation, BoundaryTriangulation}
    Γf::Union{AppendedTriangulation, CompositeTriangulation}
    E::Union{SkeletonTriangulation, Nothing}
    lateral::NamedTuple
end

struct Measures
    dΩ::Measure
    dΓ::Measure
    dΓf::Measure
    dE::Union{Measure, Nothing}
    lateral::NamedTuple
end

struct Normals
    nΓ::CellField
    nE::Union{SkeletonPair{<:CellField}, Nothing}
    lateral::NamedTuple
end

function _setup_physical_tris(cutgeo::EmbeddedDiscretization, cutgeo_facets::EmbeddedFacetDiscretization, degree::Int64)
    Ω     = Interior(cutgeo, PHYSICAL_OUT)
    Ω_act = Interior(cutgeo, ACTIVE_OUT)
    Γ     = EmbeddedBoundary(cutgeo)
    nΓ    = -get_normal_vector(Γ)
    Γf    = BoundaryTriangulation(cutgeo_facets, PHYSICAL_OUT, tags=["surface"])
    dΩ    = Measure(Ω,  degree)
    dΓ    = Measure(Γ,  degree)
    dΓf   = Measure(Γf, degree)
    return Ω, Ω_act, Γ, nΓ, Γf, dΩ, dΓ, dΓf
end

function _lateral_boundaries(model::DiscreteModel, degree::Int64, ::WallWall)
    Γw  = BoundaryTriangulation(model, tags=["walls"])
    nΓw = get_normal_vector(Γw)
    dΓw = Measure(Γw, degree)
    return (wall=Γw,), (wall=nΓw,), (wall=dΓw,)
end

function _lateral_boundaries(model::DiscreteModel, degree::Int64, ::SymmetryInlet)
    Γs  = BoundaryTriangulation(model, tags=["symmetry"])
    Γw  = BoundaryTriangulation(model, tags=["walls"])
    nΓs = get_normal_vector(Γs)
    nΓw = get_normal_vector(Γw)
    dΓs = Measure(Γs, degree)
    dΓw = Measure(Γw, degree)
    return (symmetry=Γs, wall=Γw), (symmetry=nΓs, wall=nΓw), (symmetry=dΓs, wall=dΓw)
end

function setup_triangulations(model::DiscreteModel, cutgeo::EmbeddedDiscretization, 
                                cutgeo_facets::EmbeddedFacetDiscretization, degree::Int64, 
                                ::AGFEM, lateral_tags::LateralTags)
    Ω, Ω_act, Γ, nΓ, Γf, dΩ, dΓ, dΓf = _setup_physical_tris(cutgeo, cutgeo_facets, degree)
    lat_tris, lat_norms, lat_meas       = _lateral_boundaries(model, degree, lateral_tags)
    return Triangulations(Ω, Ω_act, Γ, Γf, nothing, lat_tris), Measures(dΩ, dΓ, dΓf, nothing, lat_meas), Normals(nΓ, nothing, lat_norms)
end

function setup_triangulations(model::DiscreteModel, cutgeo::EmbeddedDiscretization, 
                                cutgeo_facets::EmbeddedFacetDiscretization, degree::Int64, 
                                ::CUTFEM, lateral_tags::LateralTags)
    Ω, Ω_act, Γ, nΓ, Γf, dΩ, dΓ, dΓf = _setup_physical_tris(cutgeo, cutgeo_facets, degree)
    E  = GhostSkeleton(cutgeo, ACTIVE_OUT)
    nE = get_normal_vector(E)
    dE = Measure(E, degree)
    lat_tris, lat_norms, lat_meas       = _lateral_boundaries(model, degree, lateral_tags)
    return Triangulations(Ω, Ω_act, Γ, Γf, E, lat_tris), Measures(dΩ, dΓ, dΓf, dE, lat_meas), Normals(nΓ, nE, lat_norms)
end

function setup_triangulations(model::DiscreteModel, cutgeo::EmbeddedDiscretization, 
                                cutgeo_facets::EmbeddedFacetDiscretization, degree::Int64, 
                                ::SBM, lateral_tags::LateralTags)
    Ω     = Interior(cutgeo, OUT)
    Γ     = Interface(Interior(cutgeo, ACTIVE_IN), Ω).⁻
    nΓ    = get_normal_vector(Γ)
    Γf    = BoundaryTriangulation(Ω, tags=["surface"])
    lat_tris, lat_norms, lat_meas = _lateral_boundaries(model, degree, lateral_tags)
    dΩ  = Measure(Ω,  degree)
    dΓ  = Measure(Γ,  degree)
    dΓf = Measure(Γf, degree)
    return Triangulations(Ω, Ω, Γ, Γf, nothing, lat_tris), Measures(dΩ, dΓ, dΓf, nothing, lat_meas), Normals(nΓ, nothing, lat_norms)
end