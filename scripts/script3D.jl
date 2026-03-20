using EmbeddedCoefficients
using Gridap
using Plots

# Physical parameters
KRs = [0.1:0.1:2.0;]
R = 9.4                         # [m]: radius
ρV = 1.0#π*R^2/2                    # [m]: area of a full horizontal cylinder (half domain)
Ks = KRs./R                     # [m⁻¹]: range of wave numbers
g = 9.81

ωs = sqrt.(Ks.*g)
n = 20

depth = 320.0       # [m]: depth


p1=plot()
p2=plot()

# Shared domain + geometry
geometry = OC3("data/meshes/oc3.stl",VectorValue(0.0,0.0,0.0))
println("domain length: ",2.0/Ks[1])
println("domain depth: ",depth)

domain   = CartesianDomain3D(2.0/Ks[1], 2.0/Ks[1], depth, (n,n,2*n))


# Run for each method and save
for method in [AGFEM(order=1)]
# for method in [AGFEM(order=1), CUTFEM(order=1), SBM(order=1), SBM(order=2)]
    added_mass = Vector{Float64}()
    added_damping = Vector{Float64}()
    params = SimulationParams(domain, geometry, method)
    setup = setup_simulation(params)
    for k in Ks
        ω = √(k * g)

        # L = max(2.0/k,50.0)
        # domain   = CartesianDomain3D(L, L, depth, (n,n,Int(trunc(depth/L)*n)))
        # println("tuple n: ",(n,n,trunc(depth/L)*n))
        # println("domain length: ",L)
        # println("domain depth: ",depth)
        # params = SimulationParams(domain, geometry, method)
        coeffs = coeff_solve(params, setup, ω, k, 1.0, g)

        display(coeffs.added_mass.*1025)
        display(coeffs.added_damping.*1025)
        # push!(added_mass,coeffs.added_mass[2,2])
        # push!(added_damping,coeffs.added_damping[2,2])
        
        # A = real(x)/(ω^2)*1025  # [kg]: added mass
        # B = imag(x)/ω*1025      # [kg/s]: added damping

        push!(added_mass,coeffs.added_mass[1,5]*1025)
        push!(added_damping,coeffs.added_damping[1,5]*1025*ω)

    end 


    plot!(p1,KRs,added_mass, label=label(method))
    plot!(p2,KRs,added_damping, label=label(method))

end
display(p1)
display(p2)