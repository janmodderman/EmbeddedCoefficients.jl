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
function coeff_solve(params::SimulationParams, ω::Float64, k::Float64, ρV::Float64)

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

    # 4. Weak forms
    a_wϕ = make_a_wϕ(params.method, k, meas, norms)
    a_wu = make_a_wu(params.method, ω, norms, meas)
    a_vϕ = make_a_vϕ(params.method, ω, norms, meas)

    # 5. Assembly + solve
    matrices = assemble_matrices(a_wϕ, a_wu, a_vϕ, spaces)
    y        = solve_system(matrices)

    # 6. Extract coefficients
    hydro_coeffs(y, ω, ρV)
end