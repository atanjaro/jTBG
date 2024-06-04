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

Returns the matrix elements of electron-phonon deformation coupling ⟨k+q|exp(iq⋅r)V(q)|k⟩.

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
function get_eph_coupling(k_el_states, kq_el_states, ephc, ph_energies, q_points, q_TF, a₁, a₂, Nk) 
    # DEBUG
    # q_points = ΓM_points
    # k_el_states = states1
    # kq_el_states = kq_states1
    # ph_energies = Ω₁_zz

    # get number of bands/orbitals
    n_bands = unit_cell.n

    # get number of branches for the phonon mode
    # energies need to be converted to matrices for this to work
    rows, n_branches = size(ph_energies)

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
        for Ω in ph_energies[:,ν]            
            prefactor = sqrt(ephc/Ω)
            push!(q_prefactors, prefactor)      
        end
    end

    
    # generate potential values for all q-points        
    for q_point in q_points
        potential = screened_coulomb_TF(q_point, q_TF, echarge)
        push!(Vqs,potential)
    end

    # get coupling matrix elements for all q-points and for all bands
    for band in 1:n_bands      
        for (r_point,q_point,V, kq_states) in zip(r_vecs, q_points, Vqs, kq_el_states)      # Nk q-points x Nk k-points
            mat_el = get_eph_matrix_elements(k_el_states, kq_states, V, q_point, r_point, band)          
            push!(matrix_elements,mat_el)       # first Nk entries are for the 1st band, second Nk entries are for the 2nd band
        end
    end

    # calculate coupling of each band to each phonon branch
    # there are either 2 or 4 branches, so this works. Maybe we can try something more generic in the future...
    if n_branches == 2
        # coupling of the first band to the first phonon branch
        for (q_point, pf, me) in zip(q_points, q_prefactors[1:Nk], matrix_elements[1:Nk])           
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        for (q_point, pf, me) in zip(q_points, q_prefactors[Nk+1:n_branches*Nk], matrix_elements[1:Nk])           
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        # coupling of the second band to the second phonon branch
        for (q_point, pf, me) in zip(q_points, q_prefactors[1:Nk], matrix_elements[Nk+1:n_branches*Nk])            
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        for (q_point, pf, me) in zip(q_points, q_prefactors[Nk+1:n_branches*Nk], matrix_elements[Nk+1:n_branches*Nk])            
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end
    elseif n_branches == 4
        # coupling of the first band to all 4 phonon branches
        for (q_point, pf, me) in zip(q_points, q_prefactors[1:Nk], matrix_elements[1:Nk])           
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        for (q_point, pf, me) in zip(q_points, q_prefactors[Nk+1:2*Nk], matrix_elements[1:Nk])            
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        for (q_point, pf, me) in zip(q_points, q_prefactors[2*Nk+1:3*Nk], matrix_elements[1:Nk])            
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        for (q_point, pf, me) in zip(q_points, q_prefactors[3*Nk+1:4*Nk], matrix_elements[1:Nk])            
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        # coupling of the second band to all 4 phonon branches
        for (q_point, pf, me) in zip(q_points, q_prefactors[1:Nk], matrix_elements[Nk+1:2*Nk])           
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        for (q_point, pf, me) in zip(q_points, q_prefactors[Nk+1:2*Nk], matrix_elements[Nk+1:2*Nk])            
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        for (q_point, pf, me) in zip(q_points, q_prefactors[2*Nk+1:3*Nk], matrix_elements[Nk+1:2*Nk])            
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end

        for (q_point, pf, me) in zip(q_points, q_prefactors[3*Nk+1:4*Nk], matrix_elements[Nk+1:2*Nk])            
            g = (pf * me * la.norm(q_point))im
            push!(g_kq, g)
        end
    end
    
    return g_kq         # for Nk q-points, [g_kq] = n_branches*Nk, with Nk couplings to Nk k-points. Each of the Nk entries per branch are coupling to each k-point i.e. the 1st entry g_kq[1] is the coupling of q₁ to k₁,k₂,...,kₙ
end




