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

    # whether information is output to terminal during runtime
    verbose = true

    # whether to plot data at the end of the simulation
    plot = true

    ##########################################
    ##          GENERAL PARAMETERS          ##
    ##########################################
    # number of k-points
    Nk = 10;

    # lattice constant (in Å)
    a = 1.0;

    # fundamental charge constant, e
    echarge = 1.0

    # reduced Planck's constant, ħ
    ħ = 1.0;

    # effective mass, m
    m_eff = 1.0;

    # Thomas-Fermi wavevector   (assuming a Thomas-Fermi screened Coulomb interaction)
    q_TF = 0.1 * m_eff * echarge / ħ^2;

    # electron-phonon coupling constants
    ephc = ħ * Nk / 2 * m_eff;      

    # initialize unit cell
    unit_cell = lu.UnitCell(
        lattice_vecs = [[3/2,√3/2]*a,
                        [3/2,-√3/2]*a],
        basis_vecs   = [[0.,0.],
                        [1.,0.]]
    )

    a₁ = unit_cell.lattice_vecs[:,1];
    a₂ = unit_cell.lattice_vecs[:,2];

    # nearest neighbor vectors
    δ₁ = [1/2,√3/2]*a;
    δ₂ = [1/2,-√3/2]*a;
    δ₃ = [-1,0]*a;
    δ_vectors = [δ₁, δ₂, δ₃];

    # electronic (Slater-Koster) parameters (in eV)
    εₚ = 0.0;
    Vppπ = -3.07;
    Vppσ = 0.0;
    slater_koster = SlaterKoster(εₚ, Vppπ, Vppσ)

    # conversion factor from cm⁻² to eV²
    conv_factor = 1/(8065.54)^2;

    # phononic parameters (in 10⁵ cm⁻² and converted to converted to eV²)
    # taken from Falkovsky, Phys. Lett. A 372 (2008) 5189–5192
    α_z = conv_factor*(-1.176e5);
    γ_z = conv_factor*(0.190e5);
    α = conv_factor*(-4.046e5);
    β = conv_factor*(1.107e5);
    γ = conv_factor*(-0.238e5);
    δ = conv_factor*(-1.096e5);
    force_constants = ForceConstants(α_z, γ_z, α, β, γ, δ)

    ########################################
    ##        BILAYER PARAMETERS          ##
    ########################################

    # TBA
    
    ###########################################
    ##           CALCULATE K-POINTS          ##
    ###########################################
    # high symmetry points
    Γ = [0,0];
    M = (unit_cell.reciprocal_vecs[:,1].+unit_cell.reciprocal_vecs[:,2])/2;
    K =  (2*unit_cell.reciprocal_vecs[:,1].+unit_cell.reciprocal_vecs[:,2])/3;                  
    K′ = (2*unit_cell.reciprocal_vecs[:,2].+unit_cell.reciprocal_vecs[:,1])/3;   

    # Generate a range of k-points along different cuts in the FBZ
    ΓM_points = [Γ + t * (M - Γ) for t in range(0, stop=1, length=Nk)];
    MK_points = [M + t * (K - M) for t in range(0, stop=1, length=Nk)];
    KΓ_points = [K + t * (Γ - K) for t in range(0, stop=1, length=Nk)];

    # # Generate (k+q)-points along cuts of the FBZ
    kq_points1 = get_all_kq_points(ΓM_points, ΓM_points);    # 10 k-points x 10 q-points. e.g. the first entry kq_points[1] corresponds to all k-points + q = (0,0)
    kq_points2 = get_all_kq_points(MK_points, MK_points);
    kq_points3 = get_all_kq_points(KΓ_points, KΓ_points);


    #####################################################
    ##          ELECTRONIC DEGREES OF FREEDOM          ##
    #####################################################      

    # generate k-dependent energies and states along all cuts of the FBZ
    (energies1, energies2, energies3, states1, states2, states3) = generate_electronic_dofs(ΓM_points, MK_points, KΓ_points, slater_koster, δ_vectors, Nk, a, false, false);

    # generate kq-dependent energies and states along all cuts of the FBZ
    (kq_energies1, kq_energies2, kq_energies3, kq_states1, kq_states2, kq_states3) = generate_electronic_dofs(kq_points1, kq_points2, kq_points3, slater_koster, δ_vectors, Nk, a, true, false);


    ###################################################
    ##          PHONONIC DEGREES OF FREEDOM          ##
    ###################################################   
    
    # generate q-dependent out-of-plane (zz) frequencies and states along all cuts of the FBZ
    (Ω₁_zz, Ω₂_zz, Ω₃_zz, zz_states1, zz_states2, zz_states3) = generate_phononic_dofs("out-of-plane", force_constants, ΓM_points, MK_points, KΓ_points, a, false);
    
    # generate q-dependent in-plane (xy) frequencies and states along all cuts of the FBZ
    (Ω₁_xy, Ω₂_xy, Ω₃_xy, xy_states1, xy_states2, xy_states3) = generate_phononic_dofs("in-plane", force_constants, ΓM_points, MK_points, KΓ_points, a, false);
    
    
    ################################################
    ##          ELECTRON-PHONON COUPLING          ##
    ################################################

    # coupling to out-of-plane modes
    coupling1 = get_eph_coupling(states1, kq_states1, ephc, Ω₁_zz, ΓM_points, q_TF, a₁, a₂, Nk);
    coupling2 = get_eph_coupling(states2, kq_states2, ephc, Ω₂_zz, MK_points, q_TF, a₁, a₂, Nk);
    coupling3 = get_eph_coupling(states3, kq_states3, ephc, Ω₃_zz, KΓ_points, q_TF, a₁, a₂, Nk);

    # coupling to in-plane modes
    coupling4 = get_eph_coupling(states1, kq_states1, ephc, Ω₁_xy, ΓM_points, q_TF, a₁, a₂, Nk);
    coupling5 = get_eph_coupling(states2, kq_states2, ephc, Ω₂_xy, MK_points, q_TF, a₁, a₂, Nk);
    coupling6 = get_eph_coupling(states3, kq_states3, ephc, Ω₃_xy, KΓ_points, q_TF, a₁, a₂, Nk);


    # collect all elements of the coupling as a function of q
    g2_kq_ΓM = []
    g2_kq_MK = []
    g2_kq_KΓ = []
    for i in 1:Nk
        ΓM_elements = [abs2.(matrix[1, 1]) for matrix in coupling1[i]][1]
        push!(g2_kq_ΓM, ΓM_elements)
        MK_elements = [abs2.(matrix[1, 1]) for matrix in coupling2[i]][1]
        push!(g2_kq_MK, MK_elements)
        KΓ_elements = [abs2.(matrix[1, 1]) for matrix in coupling3[i]][1]
        push!(g2_kq_KΓ, KΓ_elements)
    end
    
    # plot |g|² as a function of q

   


    ##########################
    ##      PLOTTING        ##
    ##########################

    if plot == true
        plot_electron_energy(energies1, energies2, energies3, ΓM_points, MK_points, KΓ_points, directory_path)

        # clear plot
        plot(fontfamily="Computer Modern",size=(800, 600),xtickfontsize=14,ytickfontsize=14,xguidefontsize=14,yguidefontsize=14)

        plot_phonon_energy(Ω₁_zz, Ω₂_zz, Ω₃_zz, Ω₁_xy, Ω₂_xy, Ω₃_xy, ΓM_points, MK_points, KΓ_points, directory_path)

    end

    return nothing
end

run()