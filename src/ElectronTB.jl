"""
    SlaterKoster( εₚ::Float64, Vppπ::Float64 )

A type defining Slater-Koster tight binding parameters.

# Fields

- `εₚ`: Carbon on-site p-orbital energy.
- `Vppπ`: hopping energy between Carbon p-orbital via π-bond.
"""
struct SlaterKoster
    # on-site p-orbital energy
    εₚ::Float64

    # p-p orbital π-bond hopping parameter
    Vppπ::Float64
end


"""
    get_exponential_sums( d_vectors::Vector{Vector{Float64}}, 
                          k_points::Vector{Vector{Float64}}, Nk::Int, a::Float64 )

Returns the sum of exponentiated dot products between k-points and distance vectors.

# Arguments

- `d_vectors::Vector{Vector{Float64}}`: Vector of nearest-neighbor distances.
- `k-points::Vector{Vector{Float64}}`: Collection of k-points in the first Brillouin Zone.
- `Nk`: Number of k-points.
- `a`: Real lattice constant.

"""
function get_exponential_sums(d_vectors, k_points, Nk, a)
    # Initialize empty arrays to store the dot products
    dot_products = [zeros(Float64, Nk) for _ in d_vectors]

    # Loop over each vector in k_points to calculate dot product with each vector in d_vectors
    for (j, d_vector) in enumerate(d_vectors)
        for i in 1:Nk
            dot_products[j][i] = la.dot(k_points[i], d_vector)
        end
    end

    # preallocate
    exp_sum = zeros(ComplexF64, Nk)

    # Loop over each k_point to calculate fk_res
    for i in 1:Nk
        exp_sum[i] = sum(exp.(a .* [dot_products[j][i] for j in 1:length(d_vectors)] .* im)) #fk(, a)
    end

    return exp_sum
end


"""
    calculate_electronic_band_structure( el_exp_sum::Vector{ComplexF64}, 
                                         εₚ::Float64, Vppπ::Float64 )

Returns the eigenenergies and eigenstates of the Graphene electronic Hamiltonian.

# Arguments

- `el_exp_sum::Vector{Float64}`: Predefined sums of exponentials.
- `slater_koster::SlaterKoster`: Slater-Koster tight binding parameters.
"""
function calculate_electronic_band_structure(el_exp_sum, slater_koster)
    energies = []
    states = []

    for i in el_exp_sum
        # Hamiltonian matrix elements
        H₁₁ = slater_koster.εₚ
        H₁₂ = slater_koster.Vppπ*i
        H₂₁ = adjoint(H₁₂)
        H₂₂ = H₁₁

        # Complete Hamiltonian matrix
        H = [H₁₁ H₁₂; H₂₁ H₂₂]

        # Diagonalize the Hamiltonian
        eigenvalues = la.eigen(H).values
        push!(energies, eigenvalues)
        
        eigenvectors = la.eigen(H).vectors
        push!(states, eigenvectors)

    end
    
    return energies, states
end

