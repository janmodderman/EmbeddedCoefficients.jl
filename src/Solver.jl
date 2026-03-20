struct SimulationParams{N}
    domain::BackgroundMesh{N}
    geometry::EmbeddedGeometry{N}
    method::EmbeddingMethod
end

struct SimulationState
    meas::Measures
    norms::Normals
    spaces::FESpaces
    dist::Union{DistanceData,Nothing}
end

function setup_simulation(params::SimulationParams)
    # 1. Domain
    model          = setup_model(params.domain)
    cut            = cut_model(model, gridap_geo(params.geometry))

    # 2. Triangulations
    degree = 2 * params.method.order
    tris, meas, norms = setup_triangulations(
        model, cut.cut, cut.facets, degree, params.method, params.domain.lateral_tag
    )
    # 3. Spaces
    spaces = setup_spaces(params.method.order, model, tris.Ω_act, cut, params.method)

    # 4. Distance functions
    dist = compute_distances(params.method, model, params.geometry, tris.Γ, degree, 0.0)
    return SimulationState(meas, norms, spaces, dist)
end     

"""
    solve(params::SimulationParams, ω::Float64, k::Float64, ρV::Float64)
        -> HydroCoefficients

Run the full pipeline for any EmbeddingMethod defined in SimulationParams.
"""
# TODO: can we pre-construct the matrices and then multiply by k, ω or products of these so we only have to assemble once?
# TODO: look into variable domain size per wavenumber k
# TODO: can we store ω in the function?
# TODO: do we want to pass ρV or calculate it inside?
function coeff_solve(params::SimulationParams, setup::SimulationState, ω::Float64, k::Float64, ρV::Float64, g::Float64)

    # 5. Weak forms
    a_wϕ = make_a_wϕ(params.method, k, ω, g, setup.meas, setup.norms, setup.dist)
    a_wu = make_a_wu(params.method, ω, params.geometry, setup.norms, setup.meas, setup.dist)
    a_vϕ = make_a_vϕ(params.method, ω, params.geometry, setup.norms, setup.meas, setup.dist)

    # 6. Assembly + solve
    matrices = assemble_matrices(a_wϕ, a_wu, a_vϕ, setup.spaces)
    y        = solve_system(matrices)

    # 7. Extract coefficients
    hydro_coeffs(y, ω, ρV)
end