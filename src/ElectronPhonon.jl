"""
    screened_coulomb_TF(q_points, q_TF, echarge,use_fft)

Given a particular q-point, generates the screened Coulomb potential matrix assuming Thomas-Fermi screening. 

"""
function screened_coulomb_TF(q_point, q_TF, echarge)

    V = 4π * echarge^2 / (la.norm(q_point)^2 + la.norm(q_TF)^2)

    return V
end


"""
get_eph_matrix_elements(k_el_states, kq_el_states, potential, q_points, r_vecs, unit_cell)

Returns the matrix elements of electron-phonon coupling ⟨k+q|exp(iq⋅r)V(q)|k⟩.

"""
function get_eph_matrix_elements(k_el_states, kq_el_states, Vq, q_point, r_point, band)
    matrix_elements = []

    for (Uk,Ukq) in zip(k_el_states, kq_el_states)
        element = transpose(Ukq[:,band]) .* Vq .* exp((la.dot(q_point,r_point))im) .* Uk[:,band]
        push!(matrix_elements, element)
    end

    return matrix_elements
end


"""
    get_eph_coupling( eph_matrix_elements::Vector{Any}, ephc::Float64, 
                    ph_energies::Matrix{AbstractFloat}, qpoints::Vector{SVector{2, Float64}} )

Computes the momentum-dependent electron-phonon coupling constant g(k,q) = i * sqrt(ħN/2MΩ(q,ν)) * ⟨k+q|exp(iq⋅r)V(q)|k⟩ * q⋅e(ν) .
"""
function get_eph_coupling(eph_matrix_elements, ephc, ph_energies, q_points)
    # get number of bands/orbitals
    n_bands = unit_cell.n

    # get number of branches for the phonon mode
    rows, n_branches = size(ph_energies1)

    # get all lattice position vectors
    r_vecs = get_positions(a₁, a₂, Nk)

    # q-dependent constants
    q_prefactors = []
    # momentum-dependent coupling
    g_kq = []
    # Coulomb potential
    Vqs = []
    # e-ph matrix elements
    matrix_elements = []
    

    # store phonon dispersion dpendendent prefactors of the form sqrt(ħN/2MΩ(q,ν))
    for ν in 1:n_branches
        for Ω in ph_energies1[:,ν]            
            prefactor = sqrt(ephc/Ω)
            push!(q_prefactors, prefactor)      # first Nk factors are for the first branch, second Nk factors are for the second branch
        end
    end

    
    # generate potential values for all q-points
    for q_point in ΓM_points
        potential = screened_coulomb_TF(q_point, q_TF, echarge)
        push!(Vqs,potential)
    end

    # get coupling matrix elements for all q-points and for all bands
    for band in 1:n_bands       
        for (r_point,q_point,V) in zip(r_vecs, ΓM_points, Vqs)      # 10 q-points x 10 k-points
            mat_el = get_eph_matrix_elements(states1, kq_states1, V, q_point, r_point, band)
            push!(matrix_elements,mat_el)       # first Nk entries are for the 1st band, second Nk entries are for the 2nd band
        end
    end

    # calculate coupling of each band to each phonon branch
    for (q_point, pf, me) in zip(ΓM_points, q_prefactors[1:10], matrix_elements)            # TODO: need to be more generic about the band number
        g = (pf * me * la.norm(q_point))im
        push!(g_kq, g)
    end

    return g_kq
end




