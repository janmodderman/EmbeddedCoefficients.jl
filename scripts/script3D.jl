using EmbeddedCoefficients
using Gridap
using Plots
using Profile

# Physical parameters
# KRs = [0.1:0.1:2.0;]
ωs = [0.2:0.2:5.0;]
g = 9.81



R = 9.4                         # [m]: radius
ρV = 1.0#π*R^2/2                    # [m]: area of a full horizontal cylinder (half domain)

KRs = ωs.^2 ./g .*R
# k = ω^2/g

Ks = KRs./R                     # [m⁻¹]: range of wave numbers
println("ωs: ",ωs)
println("Ks: ", Ks)
# ωs = sqrt.(Ks.*g)
n = 20

depth = 320.0       # [m]: depth


p1=plot()
p2=plot()

# Shared domain + geometry
geometry = OC3("data/meshes/oc3.stl",VectorValue(0.0,0.0,0.0))
println("domain length: ",1.0/Ks[1])
println("domain depth: ",depth)

domain   = CartesianDomain3D(1.0/Ks[1], 1.0/Ks[1], depth, (3*n,3*n,2*n))


# Run for each method and save
for method in [AGFEM(order=1)]
# for method in [AGFEM(order=1), CUTFEM(order=1), SBM(order=1), SBM(order=2)]
    added_mass      = Vector{Float64}()
    added_damping   = Vector{Float64}()
    params          = SimulationParams(domain, geometry, method)
    setup           = setup_simulation(params)
    matrices        = pre_assemble(params,setup)
    cache           = setup_cache(matrices)

    for k in Ks
        ω = √(k * g)
        coeffs = coeff_solve(cache, matrices, ω, k, 1.0)

        display(coeffs.added_mass.*1025)
        display(coeffs.added_damping.*1025)
        # push!(added_mass,coeffs.added_mass[2,2])
        # push!(added_damping,coeffs.added_damping[2,2])
        
        # A = real(x)/(ω^2)*1025  # [kg]: added mass
        # B = imag(x)/ω*1025      # [kg/s]: added damping

        push!(added_mass,coeffs.added_mass[1,5]*1025)
        push!(added_damping,coeffs.added_damping[1,5]*1025*ω)
        
    end 

    plot!(p1,ωs,added_mass, label=label(method))
    plot!(p2,ωs,added_damping, label=label(method))

end
display(p1)
display(p2)