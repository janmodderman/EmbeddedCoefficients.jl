using EmbeddedCoefficients
using Gridap
using GridapGmsh    
using Plots
using JSON

function conduct_study_621(HR)
    # Physical parameters
    KRs = [0.1:0.1:2.0;]
    R = 0.1                         # [m]: radius
    ρV = 2*R^2                        # [m]: area of a rectangle (half domain)
    Ks = KRs./R                     # [m⁻¹]: range of wave numbers
    g = 9.81

    p1 = plot(xlabel="k̄ [-]", ylabel="Ā₃₃ [-]")
    p2 = plot(xlabel="k̄ [-]", ylabel="B̄₃₃ [-]")

    # Shared domain + geometry
    pmid = VectorValue(0.0, 0.0)
    geometry = Rectangle(VectorValue(0.0, 0.0), 2*R, R)
    domain = GmshDomain2D("data/meshes/background_shapes.msh", lateral_tag=SymmetryInlet())

    # Dictionary to store all results for this HR value
    hr_results = Dict(
        "HR" => HR,
        "submergence_depth" => HR * R,
        "breadth" => 2*R,
        "draft" => R,
        "cylinder_center" => Dict(
            "x" => pmid[1],
            "y" => pmid[2]
        ),
        "methods" => Dict()
    )

    # Methods to study
    methods = [AGFEM(order=1), CUTFEM(order=1, γg=0.1, h=smallest_cell_size(setup_model(domain))), SBM(order=1)]

    # Run for each method and save
    for method in methods
        added_mass = Vector{Float64}()
        added_damping = Vector{Float64}()
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
            push!(added_mass, coeffs.added_mass[2, 2])
            push!(added_damping, coeffs.added_damping[2, 2])
        end

        # Store results for this method under this HR value
        hr_results["methods"][result_key] = Dict(
            "method" => method_name,
            "order" => method_order,
            "KRs" => KRs,
            "wave_numbers" => Ks,
            "added_mass" => added_mass,
            "added_damping" => added_damping,
            "metadata" => Dict(
                "num_k_values" => length(Ks),
                "k_range" => Dict("min" => minimum(Ks), "max" => maximum(Ks)),
                "added_mass_range" => Dict("min" => minimum(added_mass), "max" => maximum(added_mass)),
                "added_damping_range" => Dict("min" => minimum(added_damping), "max" => maximum(added_damping))
            )
        )

        plot!(p1, KRs, added_mass, label="$(method_name) (HR=$HR)")
        plot!(p2, KRs, added_damping, label="$(method_name) (HR=$HR)")
    end

    display(p1)
    display(p2)

    return hr_results
end

# Main execution function
function run_hr_study()
    HRs = [0.0]  # ratios of H/R as defined in the dissertation

    # Master dictionary to store all results across all HR values
    all_results = Dict(
        "metadata" => Dict(
            "description" => "Study of added mass and damping coefficients for submerged cylinder at various submergence ratios",
            "global_parameters" => Dict(
                "R" => 0.1,
                "ρV" => 0.1^2,
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
        hr_results = conduct_study_621(HR)
        
        # Store this HR's results in the master dictionary
        all_results["hr_studies"]["HR_$(HR)"] = hr_results
        
        println("Completed study for HR = $HR. Results saved.\n")
    end

    # Save complete results to JSON file
    output_file = "data/results_study_621_all_hr.json"
    mkpath(dirname(output_file))  # Create directory if it doesn't exist
    open(output_file, "w") do f
        JSON.print(f, all_results, 4)  # Pretty print with 4-space indentation
    end
    @info "Complete results saved to $output_file"

    return all_results
end

# reproduce results
run_hr_study()