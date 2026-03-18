struct Triangulations
    Ω
    Ω_act
    Γ
    Γf
    E
    lateral       
end

struct Measures
    dΩ
    dΓ
    dΓf
    dE
    lateral       
end

struct Normals
    nΓ       
    nE       
    lateral  
end

function _setup_physical_tris(cutgeo, cutgeo_facets, degree)
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

function _lateral_boundaries(model, degree, ::WallWall)
    Γw  = BoundaryTriangulation(model, tags=["walls"])
    nΓw = get_normal_vector(Γw)
    dΓw = Measure(Γw, degree)
    return (wall=Γw,), (wall=nΓw,), (wall=dΓw,)
end

function _lateral_boundaries(model, degree, ::SymmetryInlet)
    Γs  = BoundaryTriangulation(model, tags=["symmetry"])
    Γw  = BoundaryTriangulation(model, tags=["walls"])
    nΓs = get_normal_vector(Γs)
    nΓw = get_normal_vector(Γw)
    dΓs = Measure(Γs, degree)
    dΓw = Measure(Γw, degree)
    return (symmetry=Γs, wall=Γw), (symmetry=nΓs, wall=nΓw), (symmetry=dΓs, wall=dΓw)
end

function setup_triangulations(model, cutgeo, cutgeo_facets, degree, ::AGFEM, lateral_tags::LateralTags)
    Ω, Ω_act, Γ, nΓ, Γf, dΩ, dΓ, dΓf = _setup_physical_tris(cutgeo, cutgeo_facets, degree)
    lat_tris, lat_norms, lat_meas       = _lateral_boundaries(model, degree, lateral_tags)
    return Triangulations(Ω, Ω_act, Γ, Γf, nothing, lat_tris), Measures(dΩ, dΓ, dΓf, nothing, lat_meas), Normals(nΓ, nothing, lat_norms)
end

function setup_triangulations(model, cutgeo, cutgeo_facets, degree, ::CUTFEM, lateral_tags::LateralTags)
    Ω, Ω_act, Γ, nΓ, Γf, dΩ, dΓ, dΓf = _setup_physical_tris(cutgeo, cutgeo_facets, degree)
    E  = GhostSkeleton(cutgeo, ACTIVE_OUT)
    nE = get_normal_vector(E)
    dE = Measure(E, degree)
    lat_tris, lat_norms, lat_meas       = _lateral_boundaries(model, degree, lateral_tags)
    return Triangulations(Ω, Ω_act, Γ, Γf, E, lat_tris), Measures(dΩ, dΓ, dΓf, dE, lat_meas), Normals(nΓ, nE, lat_norms)
end

function setup_triangulations(model, cutgeo, cutgeo_facets, degree, ::SBM, lateral_tags::LateralTags)
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