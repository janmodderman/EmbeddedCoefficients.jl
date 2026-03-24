struct HydroCoefficients
    added_mass::Matrix{Float64}     # Mₐ
    added_damping::Matrix{Float64}  # Cₐ
end

function hydro_coeffs(y::Matrix{<:ComplexF64}, ω::Float64, ρV::Float64)
    Mₐ = real(y) / ρV / ω^2
    Cₐ = imag(y) / ρV / ω^2
    return HydroCoefficients(Mₐ, Cₐ)
end