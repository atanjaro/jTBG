import LatticeUtilities as lu
import LinearAlgebra as la
using LaTeXStrings
using Plots
using DelimitedFiles
import Base.Filesystem.mkdir
using SparseArrays 
using GeometryBasics


include("ElectronTB.jl")
include("PhononTB.jl")
include("ElectronPhonon.jl")
include("Plot.jl")
include("Utilities.jl")


function run()
    
    # specify path where files will be output
    directory_path = "/Users/andy/Documents/Codes/jTBG/output/"

    # Create the directory
    mkdir(directory_path)

    ######################################
    ##      GENERAL PARAMETERS          ##
    ######################################
    # number of k-points
    Nk = 10

    # lattice constant (in Å)
    a = 1.0

    # fundamental charge constant, e
    echarge = 1.0

    # reduced Planck's constant, ħ
    hbar = 1.0

    # effective mass, m
    m_eff = 1.0

    # Thomas-Fermi wavevector
    q_TF = 0.1 * m_eff * echarge / hbar^2

    # electron-phonon coupling constants
    ephc = hbar*Nk/2*m_eff


    # initialize unit cell
    unit_cell = lu.UnitCell(
        lattice_vecs = [[3/2,√3/2]*a,
                        [3/2,-√3/2]*a],
        basis_vecs   = [[0.,0.],
                        [1.,0.]]
    )

    a₁ = unit_cell.lattice_vecs[:,1]
    a₂ = unit_cell.lattice_vecs[:,2]

    # nearest neighbor vectors
    δ₁ = [1/2,√3/2]*a
    δ₂ = [1/2,-√3/2]*a
    δ₃ = [-1,0]*a
    δ_vectors = [δ₁, δ₂, δ₃]


    # electronic (Slater-Koster) parameters (in eV)
    εₚ = 0.0
    Vppπ = -3.07
    slater_koster = SlaterKoster(εₚ, Vppπ)

    # conversion factor from cm⁻² to eV²
    conv_factor = 1/(8065.54)^2

    # phononic parameters (in 10⁵ cm⁻² and converted to converted to eV²)
    # taken from Falkovsky, Phys. Lett. A 372 (2008) 5189–5192
    α_z = conv_factor*(-1.176e5)
    γ_z = conv_factor*(0.190e5)
    α = conv_factor*(-4.046e5)
    β = conv_factor*(1.107e5)
    γ = conv_factor*(-0.238e5)
    δ = conv_factor*(-1.096e5)
    force_constants = ForceConstants(α_z, γ_z, α, β, γ, δ)

    ########################################
    ##        BILAYER PARAMETERS          ##
    ########################################

    # (commensurate) lattice integers
    m = 9
    n = m - 1 

    # super unit cell lattice vectors
    R₁ = [-n/(m^2 + m*n + n^2), (2*m + n)/((m^2 + m*n + n^2)*√3)]
    R₂ = [m+n/(m^2 + m*n + n^2), -m+n/((m^2 + m*n + n^2)*√3)]

    # super unit cell reciprocal lattice vectors
    G₁ = n*a₁ + m*a₂
    G₂ = -m*a₁ + (n+m)*a₂

    # twist angle (in radians), set θᵣ to 0 if no twist 
    θᵣ = angle_between_vectors(m*a₁ + n*a₂,n*a₁ + m*a₂)

    # rotation matrix (-θ required to get the proper rotations!)
    R = rotation_matrix(-θᵣ)

    # generate super cell nearest neighbor vectors
    δA = [0,0]
    δB₁ = (a₁ + a₂)/3
    δB₂ = δB₁ - a₁
    δB₃ = δB₁ - a₂

    Ta₁ = R * a₁
    Ta₂ = R * a₂

    TδA = 2*(Ta₁ + Ta₂)/3
    TδB = [0,0]

    Tδ₁ = (Ta₁ + Ta₂)/3
    Tδ₂ = Tδ₁ - Ta₁
    Tδ₃ = Tδ₁ - Ta₂


    ###########################################
    ##           CALCULATE K-POINTS          ##
    ###########################################
    # high symmetry points (monolayer)
    Γ = [0,0]
    M = (unit_cell.reciprocal_vecs[:,1].+unit_cell.reciprocal_vecs[:,2])/2   
    K =  (2*unit_cell.reciprocal_vecs[:,1].+unit_cell.reciprocal_vecs[:,2])/3                  
    K′ = (2*unit_cell.reciprocal_vecs[:,2].+unit_cell.reciprocal_vecs[:,1])/3   

    # high symmetry points (bilayer)
    Mb = pi*R₁
    Kb = 2*pi ./(3*(2*R₁ + R₂))
    Kb′ = rotate_vector(Kb, 2*pi/6)

    # Generate a range of k-points along different cuts in the FBZ
    ΓM_points = [Γ + t * (M - Γ) for t in range(0, stop=1, length=Nk)]
    MK_points = [M + t * (K - M) for t in range(0, stop=1, length=Nk)]
    KΓ_points = [K + t * (Γ - K) for t in range(0, stop=1, length=Nk)]

    # # Generate (k+q)-points along cuts of the FBZ
    kq_points1 = get_all_kq_points(ΓM_points, ΓM_points)    # 10 k-points x 10 q-points. e.g. the first entry kq_points[1] corresponds to all k-points + q = (0,0)
    kq_points2 = get_all_kq_points(MK_points, MK_points)
    kq_points3 = get_all_kq_points(KΓ_points, KΓ_points)


    ##########################################################################
    ##      CALCULATE QUANTITIES FOR ELECTRONIC DEGREES OF FREEDOM          ##
    ##########################################################################       

    # generate k-dependent energies and states along all cuts of the FBZ
    (energies1, energies2, energies3, states1, states2, states3) = generate_electronic_dofs(ΓM_points, MK_points, KΓ_points, slater_koster, δ_vectors, Nk, a, false, false)

    # generate kq-dependent energies and states along all cuts of the FBZ
    (kq_energies1, kq_energies2, kq_energies3, kq_states1, kq_states2, kq_states3) = generate_electronic_dofs(kq_points1, kq_points2, kq_points3, slater_koster, δ_vectors, Nk, a, true, false)


    ########################################################################
    ##      CALCULATE QUANTITIES FOR PHONONIC DEGREES OF FREEDOM          ##
    ########################################################################    
    # get out-of-plane phonon dispersion
    (ph_energies1,phstates1) = calculate_out_of_plane_phononic_dispersion(force_constants, ΓM_points,a)
    (ph_energies2,phstates2) = calculate_out_of_plane_phononic_dispersion(force_constants, MK_points,a)
    (ph_energies3,phstates3) = calculate_out_of_plane_phononic_dispersion(force_constants, KΓ_points,a)

    # get in-plane phonon dipsersion
    (ph_energies4,phstates4) = calculate_in_plane_phononic_dispersion(force_constants, ΓM_points,a)
    (ph_energies5,phstates5) = calculate_in_plane_phononic_dispersion(force_constants, MK_points,a)
    (ph_energies6,phstates6) = calculate_in_plane_phononic_dispersion(force_constants, KΓ_points,a)

    # convert phoninic energies to matrices
    ph_energies1 = permutedims(hcat(ph_energies1...))
    ph_energies2 = permutedims(hcat(ph_energies2...))
    ph_energies3 = permutedims(hcat(ph_energies3...))
    ph_energies4 = permutedims(hcat(ph_energies4...))
    ph_energies5 = permutedims(hcat(ph_energies5...))
    ph_energies6 = permutedims(hcat(ph_energies6...))

    # write energies to file
    writedlm(directory_path*"gamma_M_zz_ph_energies.csv", ph_energies1)
    writedlm(directory_path*"M_K_zz_ph_energies.csv", ph_energies2)
    writedlm(directory_path*"K_gamma_zz_ph_energies.csv", ph_energies3)
    writedlm(directory_path*"gamma_M_xy_ph_energies.csv", ph_energies4)
    writedlm(directory_path*"M_K_zz_xy_energies.csv", ph_energies5)
    writedlm(directory_path*"K_gamma_xy_ph_energies.csv", ph_energies6)

    # convert phononic states to matrices
    phstates1 = permutedims(hcat(phstates1...))
    phstates2 = permutedims(hcat(phstates2...))
    phstates3 = permutedims(hcat(phstates3...))
    phstates4 = permutedims(hcat(phstates4...))
    phstates5 = permutedims(hcat(phstates5...))
    phstates6 = permutedims(hcat(phstates6...))

    # write states to file
    writedlm(directory_path*"gamma_M_zz_ph_states.csv", phstates1)
    writedlm(directory_path*"M_K_zz_ph_states.csv", phstates2)
    writedlm(directory_path*"K_gamma_zz_ph_states.csv", phstates3)
    writedlm(directory_path*"gamma_M_xy_ph_states.csv", phstates4)
    writedlm(directory_path*"M_K_zz_xy_states.csv", phstates5)
    writedlm(directory_path*"K_gamma_xy_ph_states.csv", phstates6)

    #####################################################################
    ##      CALCULATE QUANTITIES FOR ELECTRON-PHONON COUPLING          ##
    #####################################################################

    # coupling to out-of-plane modes, (q = 0) for now
    coupling1 = get_eph_coupling(ephc, ph_energies1, ΓM_points, q_TF, a₁, a₂, Nk)
    coupling2 = get_eph_coupling(ephc, ph_energies1, MK_points, q_TF, a₁, a₂, Nk)
    coupling3 = get_eph_coupling(ephc, ph_energies1, KΓ_points, q_TF, a₁, a₂, Nk)

    # collect all diagonal elements of the coupling
    AA_elements = [matrix[1, 1] for matrix in coupling1[2]]
    BB_elements = [matrix[2, 2] for matrix in coupling1[2]]

    # Display the result
    println(AA_elements)



    ##########################
    ##      PLOTTING        ##
    ##########################
    
    plot_electron_energy(el_energies1,el_energies2,el_energies3,ΓM_points,MK_points,KΓ_points,directory_path)

    # clear plot
    plot(fontfamily="Computer Modern",size=(800, 600),xtickfontsize=14,ytickfontsize=14,xguidefontsize=14,yguidefontsize=14)

    plot_phonon_energy(ph_energies1,ph_energies2,ph_energies3,ph_energies4,ph_energies5,ph_energies6,ΓM_points,MK_points,KΓ_points,directory_path)

    return nothing
end

run()