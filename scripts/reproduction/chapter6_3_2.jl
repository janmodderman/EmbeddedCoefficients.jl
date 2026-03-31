using EmbeddedCoefficients
using Gridap
using GridapGmsh    
using Plots
using JSON

function conduct_study_632(HR)
    # Physical parameters
    KRs = [0.1:0.1:2.0;]
    R = 0.1                         # [m]: radius
    ρV = 1.0                    # [m]: don't normalize
    Ks = KRs./R                     # [m⁻¹]: range of wave numbers
    g = 9.81

    # Define the motion directions (matrix indices) to track
    # Format: (row, col) => "description"
    motion_directions = Dict(
        (1, 1) => "surge",
        (2, 2) => "sway",
        (3, 3) => "heave",
        (4, 4) => "roll",
        (5, 5) => "pitch",
        (6, 6) => "yaw",
        (2, 4) => "sway-roll",
        (4, 2) => "roll-sway",
        (1, 5) => "surge-pitch",
        (5, 1) => "pitch-surge"
    )

    # Initialize storage for all motion directions
    added_mass_data = Dict(dir => Vector{Float64}() for dir in keys(motion_directions))
    added_damping_data = Dict(dir => Vector{Float64}() for dir in keys(motion_directions))

    # Create plots for each direction
    plots = Dict()
    for (dir, description) in motion_directions
        p_mass = plot(xlabel="k̄ [-]", ylabel="A$(dir[1])$(dir[2]) [-]", title="Added Mass: $description")
        p_damp = plot(xlabel="k̄ [-]", ylabel="B$(dir[1])$(dir[2]) [-]", title="Added Damping: $description")
        plots[dir] = (mass=p_mass, damping=p_damp)
    end

    # Shared domain + geometry
    pmid = VectorValue(0.0, 0.0, 0.0)
    geometry = OC3("data/meshes/oc3.stl", pmid)
    domain = GmshDomain3D("data/meshes/background3doc3_circ.msh", lateral_tag=WallWall())

    # Dictionary to store all results for this HR value
    hr_results = Dict(
        "HR" => HR,
        "submergence_depth" => 120.0,
        "upper radius" => 3.25,
        "lower radius" => 4.70,
        "cylinder_center" => Dict(
            "x" => pmid[1],
            "y" => pmid[2],
            "z" => pmid[3]
        ),
        "motion_directions" => motion_directions,
        "methods" => Dict()
    )

    # Methods to study
    methods = [AGFEM(order=1)]

    # Run for each method and save
    for method in methods
        # Initialize storage for this method
        method_added_mass = Dict(dir => Vector{Float64}() for dir in keys(motion_directions))
        method_added_damping = Dict(dir => Vector{Float64}() for dir in keys(motion_directions))

        params = SimulationParams(domain, geometry, method)
        setup = setup_simulation(params)
        matrices = pre_assemble(params, setup)
        cache = setup_cache(matrices)

        # Get method name and order for organizing results
        method_name = label(method)  # e.g., "AGFEM", "CUTFEM", "SBM"
        method_order = method.order   # Get the order (assumes all methods have .order field)
        result_key = "$(method_name)_order_$(method_order)"

        for k in Ks
            ω = √(k * g)
            coeffs = coeff_solve(cache, matrices, ω, k, ρV)
            
            # Extract all motion directions
            for (i, j) in keys(motion_directions)
                push!(method_added_mass[(i, j)], coeffs.added_mass[i, j])
                push!(method_added_damping[(i, j)], coeffs.added_damping[i, j])
            end
        end

        # Store results for this method under this HR value
        hr_results["methods"][result_key] = Dict(
            "method" => method_name,
            "order" => method_order,
            "KRs" => KRs,
            "wave_numbers" => Ks,
            "motion_directions" => Dict(
                "added_mass" => Dict(
                    string(dir) => values for (dir, values) in method_added_mass
                ),
                "added_damping" => Dict(
                    string(dir) => values for (dir, values) in method_added_damping
                )
            ),
            "metadata" => Dict(
                "num_k_values" => length(Ks),
                "k_range" => Dict("min" => minimum(Ks), "max" => maximum(Ks)),
                "added_mass_ranges" => Dict(
                    string(dir) => Dict("min" => minimum(values), "max" => maximum(values))
                    for (dir, values) in method_added_mass
                ),
                "added_damping_ranges" => Dict(
                    string(dir) => Dict("min" => minimum(values), "max" => maximum(values))
                    for (dir, values) in method_added_damping
                )
            )
        )

        # Plot results for each motion direction
        for (dir, description) in motion_directions
            plot!(plots[dir].mass, KRs, method_added_mass[dir], label="$(method_name) (HR=$HR)")
            plot!(plots[dir].damping, KRs, method_added_damping[dir], label="$(method_name) (HR=$HR)")
        end
    end

    # Display plots for each motion direction
    for (dir, description) in motion_directions
        display(plots[dir].mass)
        display(plots[dir].damping)
    end

    return hr_results
end

# Main execution function
function run_hr_study()
    HRs = [0.0]  # ratios of H/R as defined in the dissertation

    # Master dictionary to store all results across all HR values
    all_results = Dict(
        "metadata" => Dict(
            "description" => "Study of added mass and damping coefficients for OC3",
            "global_parameters" => Dict(
                "R" => 0.1,
                "ρV" => 0.1,
                "g" => 9.81,
                "KRs" => collect(0.1:0.1:2.0),
                "num_wavelengths" => 20
            ),
            "HRs_studied" => HRs,
            "num_hr_values" => length(HRs)
        ),
        "hr_studies" => Dict()
    )

    # Run study for each HR value
    for HR in HRs
        println("Running study for HR = $HR...")
        hr_results = conduct_study_632(HR)
        
        # Store this HR's results in the master dictionary
        all_results["hr_studies"]["HR_$(HR)"] = hr_results
        
        println("Completed study for HR = $HR. Results saved.\n")
    end

    # Save complete results to JSON file
    output_file = "data/results_study_632_all_hr.json"
    mkpath(dirname(output_file))  # Create directory if it doesn't exist
    open(output_file, "w") do f
        JSON.print(f, all_results, 4)  # Pretty print with 4-space indentation
    end
    @info "Complete results saved to $output_file"

    return all_results
end

# reproduce results
run_hr_study()