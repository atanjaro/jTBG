import LatticeUtilities as lu
import LinearAlgebra as la

include("../src/ElectronTB.jl")
include("../src/PhononTB.jl")
include("../src/ElectronPhonon.jl")
include("../src/Utilities.jl")

# check calculation of graphene dispersion
Nk = 10
a = 1.0
unit_cell = lu.UnitCell(
        lattice_vecs = [[3/2,√3/2]*a,
                        [3/2,-√3/2]*a],
        basis_vecs   = [[0.,0.],
                        [1.,0.]]
    )
δ₁ = [1/2,√3/2]*a
δ₂ = [1/2,-√3/2]*a
δ₃ = [-1,0]*a
δ_vectors = [δ₁, δ₂, δ₃]
εₚ = 0.0
Vppπ = -3.07
slater_koster = SlaterKoster(εₚ, Vppπ)

Γ = [0,0]
M = (unit_cell.reciprocal_vecs[:,1].+unit_cell.reciprocal_vecs[:,2])/2  
ΓM_points = [Γ + t * (M - Γ) for t in range(0, stop=1, length=Nk)] 

el_exp_sum = get_exponential_sums(δ_vectors, ΓM_points, Nk, a)
(el_energies, states) = calculate_electronic_band_structure(el_exp_sum1, slater_koster)
