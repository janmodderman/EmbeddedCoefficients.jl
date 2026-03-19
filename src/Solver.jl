struct SimulationParams{N}
    domain::BackgroundMesh{N}
    geometry::EmbeddedGeometry{N}
    method::EmbeddingMethod
end

"""
    solve(params::SimulationParams, ω::Float64, k::Float64, ρV::Float64)
        -> HydroCoefficients

Run the full pipeline for any EmbeddingMethod defined in SimulationParams.
"""
# TODO: pull out domain + triangulations + spaces + compute_distances so we only have to run it once? 
# TODO: look into variable domain size per wavenumber k
# TODO: can we store ω in the function?
# TODO: do we want to pass ρV or calculate it inside?
function coeff_solve(params::SimulationParams, ω::Float64, k::Float64, ρV::Float64, g::Float64)

    # 1. Domain
    model          = setup_model(params.domain)
    cutgeo, facets = cut_model(model, gridap_geo(params.geometry))

    # 2. Triangulations
    degree = 2 * params.method.order
    tris, meas, norms = setup_triangulations(
        model, cutgeo, facets, degree, params.method, params.domain.lateral_tag
    )
    # 3. Spaces
    spaces = setup_spaces(params.method.order, model, tris.Ω_act, cutgeo, params.method)

    # 4. Distance functions
    dist = compute_distances(params.method, model, params.geometry, tris.Γ, degree, 0.0)

    # 5. Weak forms
    a_wϕ = make_a_wϕ(params.method, k, ω, g, meas, norms, dist)
    a_wu = make_a_wu(params.method, ω, params.geometry, norms, meas, dist)
    a_vϕ = make_a_vϕ(params.method, ω, params.geometry, norms, meas, dist)

    # 6. Assembly + solve
    matrices = assemble_matrices(a_wϕ, a_wu, a_vϕ, spaces)
    y        = solve_system(matrices)

    # 7. Extract coefficients
    hydro_coeffs(y, ω, ρV)
end