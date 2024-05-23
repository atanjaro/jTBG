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


"""
    generate_electronic_dofs(k_points1, k_points2, k_points3, slater_koster, δ_vectors, Nk, a, directory_path, k_name)

Calculates the electronic band structure and writes quantities to file. Returns the energies
and associated states.

"""
function generate_electronic_dofs(k_points1, k_points2, k_points3, slater_koster, δ_vectors, Nk, a, kq_flag, write)
    if kq_flag == true
        # pre-allocate sums of exponentials
        exp_sums1 = []
        exp_sums2 = []
        exp_sums3 = []
        # generate an array of k+q exponntial sums. For instance, the first entry of exp_sums corresponds to the sums assocaited with k + (0,0)
        for i in 1:Nk
            exp_sum1 = get_exponential_sums(δ_vectors, k_points1[i], Nk, a)  
            exp_sum2 = get_exponential_sums(δ_vectors, k_points2[i], Nk, a)
            exp_sum3 = get_exponential_sums(δ_vectors, k_points3[i], Nk, a)
            push!(exp_sums1, exp_sum1)
            push!(exp_sums2, exp_sum2)
            push!(exp_sums3, exp_sum3)
        end

        all_energies1 = []
        all_energies2 = []
        all_energies3 = []
        all_states1 = []
        all_states2 = []
        all_states3 = []
        for i in 1:Nk
            # get (k+q) band structure 
            (energies1, states1) = calculate_electronic_band_structure(exp_sums1[i], slater_koster)
            (energies2, states2) = calculate_electronic_band_structure(exp_sums2[i], slater_koster)
            (energies3, states3) = calculate_electronic_band_structure(exp_sums3[i], slater_koster)
            push!(all_energies1, energies1)
            push!(all_energies2, energies2)
            push!(all_energies3, energies3)
            push!(all_states1, states1)
            push!(all_states2, states2)
            push!(all_states3, states3)
        end

        if write == true
            # write energies to file
            writedlm(directory_path*"electronic_kq_energies1.csv", all_energies1)
            writedlm(directory_path*"electronic_kq_energies2.csv", all_energies2)
            writedlm(directory_path*"electronic_kq_energies3.csv", all_energies3)

            # write states to file
            writedlm(directory_path*"electronic_kq_states1.csv", all_states1)
            writedlm(directory_path*"electronic_kq_states2.csv", all_states2)
            writedlm(directory_path*"electronic_kq_states3.csv", all_states3)
        end

        return all_energies1, all_energies2, all_energies3, all_states1, all_states2, all_states3

    elseif kq_flag == false
        # pre-allocate sums of exponentials  
        exp_sum1 = get_exponential_sums(δ_vectors, k_points1, Nk, a)
        exp_sum2 = get_exponential_sums(δ_vectors, k_points2, Nk, a)
        exp_sum3 = get_exponential_sums(δ_vectors, k_points3, Nk, a)

        # get band structure
        (energies1, states1) = calculate_electronic_band_structure(exp_sum1, slater_koster)
        (energies2, states2) = calculate_electronic_band_structure(exp_sum2, slater_koster)
        (energies3, states3) = calculate_electronic_band_structure(exp_sum3, slater_koster)

        # convert electronic energies to matrices 
        energies1 = permutedims(hcat(energies1...))
        energies2 = permutedims(hcat(energies2...))
        energies3 = permutedims(hcat(energies3...))

        if write == true
            # write energies to file
            writedlm(directory_path*"electronic_k_energies1.csv", energies1)
            writedlm(directory_path*"electronic_k_energies2.csv", energies2)
            writedlm(directory_path*"electronic_k_energies3.csv", energies3)

            # write states to file
            writedlm(directory_path*"electronic_k_states1.csv", states1)
            writedlm(directory_path*"electronic_k_states2.csv", states2)
            writedlm(directory_path*"electronic_k_states3.csv", states3)
        end

        return energies1, energies2, energies3, states1, states2, states3
    end
end


"""
    generate_electronic_dofs(slater_koster, δ_vectors, ΓM_points, Nk, a)

Calculates the electronic band structure and writes quantities to file. Returns the energies
and associated states.

"""
function generate_electronic_dofs(kq_points1, kq_points2, kq_points3, slater_koster, δ_vectors, Nk, a, directory_path, kq_name)


    # get (k+q) band structure 
    (energies1, kq_states1) = calculate_electronic_band_structure(kq_el_exp_sum1, slater_koster)
    (kq_el_energies2, kq_states2) = calculate_electronic_band_structure(kq_el_exp_sum2, slater_koster)
    (kq_el_energies3, kq_states3) = calculate_electronic_band_structure(kq_el_exp_sum3, slater_koster)

    # convert (k+q) energies to matrices 
    kq_el_energies1 = permutedims(hcat(kq_el_energies1...))
    kq_el_energies2 = permutedims(hcat(kq_el_energies2...))
    kq_el_energies3 = permutedims(hcat(kq_el_energies3...))

    # write energies to file
    writedlm(directory_path*k_name*"_electronic_energies1.csv", energies1)
    writedlm(directory_path*k_name*"_electronic_energies2.csv", energies2)
    writedlm(directory_path*k_name*"_electronic_energies3.csv", energies3)

    # write states to file
    writedlm(directory_path*k_name*"_electronic_states1.csv", states1)
    writedlm(directory_path*k_name*"_electronic_states2.csv", states2)
    writedlm(directory_path*k_name*"_electronic_states3.csv", states3)





    return energies1, energies2, energies3, states1, states2, states3
end


