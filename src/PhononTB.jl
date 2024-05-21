"""
    ForceConstants( α_z::Float64, γ_z::Float64, α::Float64, 
                    β::Float64, γ::Float64, δ::Float64 )

A type defining dynamical matrix force constants.

"""
struct ForceConstants
    # out-of-plane force constants
    α_z::Float64
    γ_z::Float64
    # in-plane force constants
    α::Float64      
    β::Float64      
    γ::Float64
    δ::Float64
end


"""
    calculate_out_of_plane_phononic_dispersion( force_constants::ForceConstants, 
                                                q_points::Vector{Vector{Float64}}, a::Float64 )

Returns the eignenergies and eigenstates of the Graphene out-of-plane (zz) phonon modes.

"""
function calculate_out_of_plane_phononic_dispersion(force_constants, q_points, a)
    energies = []
    states = []


    for q in q_points
        # dynamical matrix elements
        D₁ = Φ₁_z(force_constants.α_z, force_constants.γ_z, q,a)                
        D₂ = Φ₂_z(force_constants.α_z, q,a)               
        D₃ = adjoint(D₂)
        D₄ = D₁

        # complete dynamical matrix
        D = [D₁ D₂; D₃ D₄]

        # Diagonalize out-of-plane dynamical matrix
        eigenvalues = la.eigen(D).values
        push!(energies, sqrt.(eigenvalues))
                    
        eigenvectors = la.eigen(D).vectors
        push!(states, eigenvectors)
    end

    return energies, states
end

"""
    calculate_in_plane_phononic_dispersion( force_constants::ForceConstants, 
                                            q_points::Vector{Vector{Float64}}, a::Flaot64 )

Returns the eignenergies and eigenstates of the Graphene in-plane (xx, xy, yx, yy) phonon modes.

"""
function calculate_in_plane_phononic_dispersion(force_constants, q_points, a)
    energies = []
    states = []

    for q in q_points
        D₁ = Φ₁(force_constants.α,force_constants.γ,q,a)
        D₂ = Φ₂(force_constants.δ,q,a)
        D₃ = Φ₃(force_constants.α,q,a)
        D₄ = Φ₄(force_constants.β,q,a)
        D₅ = Φ₅(force_constants.β,q,a)
        D₆ = Φ₆(force_constants.δ,q,a)

        D = [D₁ D₂ D₃ D₄; adjoint(D₂) D₁ D₅ D₃; adjoint(D₃) adjoint(D₅) D₁ D₆; adjoint(D₄) adjoint(D₃) adjoint(D₆) D₁]

        # Diagonalize out-of-plane dynamical matrix
        eigenvalues = la.eigen(D).values
        push!(energies, sqrt.(eigenvalues))
                    
        eigenvectors = la.eigen(D).vectors
        push!(states, eigenvectors)
    end
    return energies, states
end


# convenience functions for evaluating out-of-plane dispersion
function Φ₁_z(α_z, γ_z, q, a)
    D₁ = 2*γ_z*(cos(a*sqrt(3)*q[2])+2*cos(3*a*q[1]/2)*cos(a*sqrt(3)*q[2]/2)-3)-3*α_z
    return D₁
end

function Φ₂_z(α_z, q, a)
    arg = q[1] 
    D₂ = α_z*(exp((a*arg)im)+2*exp(-(a*arg/2)im)*cos(sqrt(3)*a*q[2]/2))
    return D₂
end

# convenience functions for evaluating in-plane dispersion
function Φ₁(α,γ,q,a) 
    arg = a* sqrt(3) * q[2]
    D₁ = 2 * γ * (cos(arg) + 2 * cos(3 * a * q[1] / 2) * cos(arg / 2) - 3) - 3 * α
    return D₁
end

function Φ₂(δ,q,a)
    arg = a* sqrt(3) * q[2]
    D₂ = δ * (exp((arg)im) + 2 * cos(3 * a * q[1] / 2 + 2 * pi / 3) * exp(-(arg / 2)im)) +
        adjoint(δ) * (exp(-(arg)im) + 2 * cos(3 * a * q[1] / 2 - 2 * pi / 3) * exp((arg / 2)im))
    return D₂
end

function Φ₃(α,q,a)
    D₃ = α * (exp((a * q[1])im) + 2 * cos(a * sqrt(3) * q[2] / 2) * exp(-(a * q[1] / 2)im))
    return D₃
end

function Φ₄(β,q,a) 
    D₄ = β * (exp((a * q[1])im) + 2 * cos(a * sqrt(3) * q[2] / 2 - 2 * pi / 3) * exp(-(a * q[1] / 2)im))
    return D₄
end

function Φ₅(β,q,a) 
    D₅ = adjoint(β) * (exp((a * q[1])im) + 2 * cos(a * sqrt(3) * q[2] / 2 + 2 * pi / 3) * exp(-(a * q[1] / 2)im))
    return D₅
end

function Φ₆(δ,q,a) 
    arg = a * sqrt(3) * q[2]
    D₆ = δ * (exp(-(arg)im) + 2 * cos(3 * a * q[1] / 2 - 2 * pi / 3) * exp((arg / 2)im)) +
        adjoint(δ) * (exp((arg)im) + 2 * cos(3 * a * q[1] / 2 + 2 * pi / 3) * exp(-(arg / 2)im))
    return D₆
end
