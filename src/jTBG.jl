module jTBG

import LinearAlgebra as la

include("ElectronTB.jl")
export SlaterKoster
export get_exponential_sums
export calculate_electronic_band_structure

include("PhononTB.jl")
export ForceConstants
export calculate_out_of_plane_phononic_dispersion
export calculate_in_plane_phononic_dispersion

end # module