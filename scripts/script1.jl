using EmbeddedCoefficients
using Gridap
using Plots

# Physical parameters
KRs = [0.1:0.1:2.0;]
R = 0.1                         # [m]: radius
ρV = π*R^2/2                    # [m]: area of a full horizontal cylinder (half domain)
Ks = KRs./R                     # [m⁻¹]: range of wave numbers
g = 9.81
n = 350

p1=plot()
p2=plot()

# Shared domain + geometry
geometry = Circle(VectorValue(0.0, 0.0), R)
println("domain length: ",1.0/KRs[1])
println("domain depth: ",20*R)
domain   = CartesianDomain2D(1.0/KRs[1], 20*R, (5*n,n), lateral_tag=SymmetryInlet())
# domain   = CartesianDomain2D(1.0/KRs[1], 20*R, (5*n,n), lateral_tag=WallWall())


# Run for each method and save
# for method in [SBM(order=1)]
for method in [AGFEM(order=1), CUTFEM(order=1), SBM(order=1), SBM(order=2)]
    added_mass = Vector{Float64}()
    added_damping = Vector{Float64}()
    for k in Ks
        ω = √(k * g)

        # domain   = CartesianDomain2D(1.0/(k/R), 20*R, (50,50), lateral_tag=SymmetryInlet())

        params = SimulationParams(domain, geometry, method)
        coeffs = coeff_solve(params, ω, k, ρV, g)

        # display(coeffs.added_mass)
        # display(coeffs.added_damping)
        push!(added_mass,coeffs.added_mass[2,2])
        push!(added_damping,coeffs.added_damping[2,2])
        

        # push!(added_mass,coeffs.added_mass[3,3])
        # push!(added_damping,coeffs.added_damping[3,3])

    # name       = string(typeof(method))
    # output_dir = joinpath(@__DIR__, "../../output/2d_circle/$(name)")
    # mkpath(output_dir)

    # df = DataFrame(
    #     added_mass    = vec(coeffs.added_mass),
    #     added_damping = vec(coeffs.added_damping)
    # )
    # CSV.write(joinpath(output_dir, "coefficients.csv"), df)

    # println("$(name) — Mₐ: $(coeffs.added_mass), Cₐ: $(coeffs.added_damping)")
    end 


    plot!(p1,KRs,added_mass, label=label(method))
    plot!(p2,KRs,added_damping, label=label(method))

end
display(p1)
display(p2)