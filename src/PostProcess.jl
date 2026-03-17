struct HydroCoefficients
    added_mass::Matrix{Float64}     # Mₐ
    added_damping::Matrix{Float64}  # Cₐ
end

"""
Extract non-dimensionalised added mass and damping from the y matrix.

# Arguments
- `y`  : output of solve_system
- `ω`  : oscillation frequency [rad/s]
- `ρV` : ρ * V (fluid density × reference volume)
"""
function hydro_coeffs(y, ω::Float64, ρV::Float64)
    Mₐ = real(y) / ρV / ω^2
    Cₐ = imag(y) / ρV / ω^2
    HydroCoefficients(Mₐ, Cₐ)
end