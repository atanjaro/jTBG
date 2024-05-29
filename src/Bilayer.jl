using Base.Iterators: product

include("ElectronTB.jl")
include("PhononTB.jl")
include("ElectronPhonon.jl")
include("Plot.jl")
include("Utilities.jl")


# tight-binding parameters
Vppπ = 2.8 # nearest neighbor hopping in eV
Vppσ = 0.0  # next-nearest neighbor hopping in eV 

# (commensurate) lattice integers
m = 9
n = m - 1 

# total number of sites in the extended (super) unit cell
N_suc = n^2 + n*m + m^2

# extended unit cell lattice vectors
R₁ = [-n/(m^2 + m*n + n^2), (2*m + n)/((m^2 + m*n + n^2)*√3)]
R₂ = [m+n/(m^2 + m*n + n^2), -m+n/((m^2 + m*n + n^2)*√3)]

# extended unit cell reciprocal lattice vectors
G₁ = n*a₁ + m*a₂
G₂ = -m*a₁ + (n+m)*a₂

# twist angle (in radians), set θᵣ to 0 if no twist 
# TODO: there is a bug which causing the incorrect angle to be reported
θᵣ = angle_between_vectors(m*a₁ + n*a₂,n*a₁ + m*a₂)  
# convert angle to degrees
θ_deg = (180*θᵣ/pi)

# rotation matrix 
R = rotation_matrix(-θᵣ)

# generate nearest neighbor vectors for the extended cell
δA = [0,0]
δB₁ = (a₁ + a₂)/3
δB₂ = δB₁ - a₁
δB₃ = δB₁ - a₂

# twisted layer extended unit cell lattice vectors
Ta₁ = R * a₁
Ta₂ = R * a₂

# nearest neighbor vectors for the twisted layer
TδA = 2*(Ta₁ + Ta₂)/3
TδB = [0,0]
Tδ₁ = (Ta₁ + Ta₂)/3
Tδ₂ = Tδ₁ - Ta₁
Tδ₃ = Tδ₁ - Ta₂

# high symmetry points (bilayer)
Mb = pi*R₁
Kb = 2*pi ./(3*(2*R₁ + R₂))
Kb′ = rotate_vector(Kb, 2*pi/6)

"""
    generate_extended_unit_cell(n, m, G₁, G₂, δB₁, δB₂, δB₃, layer)

Generates list of sites and bonds in the extended unit cell of either the untwisted
or untwisted layer.

"""
function generate_extended_unit_cell(n, m, G₁, G₂, δB₁, δB₂, δB₃, layer)
    # where points and bonds will be stored
    A_points = Vector{Any}()
    B_points = Vector{Any}()
    bonds = Vector{Any}()

    # generate all sites in the extended unit cell of the untwisted layer
    if layer == 1
        # A points and bonds
        for (i, j) in product(-2*(n + m):2*(n + m), -2*(n + m):2*(n + m))
            if region_member_parallelogram([0, 0], G₁ * 0.9999999, G₂ * 0.9999999, rA(i, j))
                push!(A_points, (i, j, rA(i, j)))
                push!(bonds, (rA(i, j), rA(i, j) + δB₁))
                push!(bonds, (rA(i, j), rA(i, j) + δB₂))
                push!(bonds, (rA(i, j), rA(i, j) + δB₃))
            end
        end

        # B points
        for (i, j) in product(-2*(n + m):2*(n + m), -2*(n + m):2*(n + m))
            if region_member_parallelogram([0, 0], G₁, G₂, rB(i, j))
                push!(B_points, (i, j, rB(i, j)))
            end
        end

        return A_points, B_points, bonds
        
    # generate all sites in the extended unit cell for the twisted layer    
    elseif layer == 2
        # A points and bonds
        for (i, j) in product(-2*(n + m):2*(n + m), -2*(n + m):2*(n + m))
            if region_member_parallelogram([0, 0], G₁ * 0.9999999, G₂ * 0.9999999, TrA(i, j))
                push!(A_points, (i, j, TrA(i, j)))
                push!(bonds, (TrA(i, j), TrA(i, j) + Tδ₁))
                push!(bonds, (TrA(i, j), TrA(i, j) + Tδ₂))
                push!(bonds, (TrA(i, j), TrA(i, j) + Tδ₃))
            end
        end

        # B points
        for (i, j) in product(-2*(n + m):2*(n + m), -2*(n + m):2*(n + m))
            if region_member_parallelogram([0, 0], G₁, G₂, rB(i, j))
                push!(B_points, (i, j, TrB(i, j)))
            end
        end
        
        return A_points, B_points, bonds
    else
        println("ERROR: NOT A VALID LAYER!")
        println("Enter either layer '1' (untwisted) or '2' (twisted)")
        return nothing
    end
end

# get all points of the untwisted layer (layer 1)
A₁_points, B₁_points, bonds = generate_extended_unit_cell(n, m, G₁, G₂, δB₁, δB₂, δB₃, 1);

# get all points of the twisted layer (layer 2)
A₂_points, B₂_points, Tbonds = generate_extended_unit_cell(n, m, G₁, G₂, δB₁, δB₂, δB₃, 2);

# collect all lattice points
lattice_points = vcat(A₁_points, B₁_points, A₂_points, B₂_points);



# TODO: while running over all of the sites for intralayer hopping, we need to account for the k-points
# that is, kx and ky at Γ and K
function generate_intralayer_hopping_hamiltonian(slater_koster, N_suc, lattice_points, k_point, G₁, G₂, δB₁, δB₂, δB₃, layer)
    k_point = Kb
    
    Vppπ = slater_koster.Vppπ
    Vppσ = slater_koster.Vppσ

    δ_to_points = [rotation_matrix(k * 2 * pi / 6) * a₁ for k in 0:2]
    Tδ_to_points = [rotation_matrix(k * 2 * pi / 6) * Ta₁ for k in 0:2]

    intralayer_connections = []

    kx = k_point[1]
    ky = k_point[2]

    if layer == 1
        # run over A sublattice of the untwisted (1st) layer
        for x₁ in 1:N_suc
            # Position of A site in 1st layer
            v₀ = lattice_points[x₁][3]

            # First nearest neighbor
            v₁ = v₀ + δB₁
            while !region_member_parallelogram([0, 0], G₁ * 0.999999, G₂ * 0.999999, v₁)
                α = rand(-1:1)
                β = rand(-1:1)
                v₁ = v₀ + δB₁ + α * G₁ + β * G₂
            end
            position_1 = findfirst(x -> x == v₁, lattice_points[1:2*N_suc])
            if position_1 !== nothing
                push!(intralayer_connections, (x₁, position_1) => -Vppπ * exp(1im * [kx, ky] .* δB₁))         
                push!(intralayer_connections, (position_1, x₁) => -Vppπ * exp(-1im * [kx, ky] .* δB₁))      
            end
        end

        for x₁ in 1:N_suc
            # Second nearest neighbor
            v₁ = v₀ + δB₂
            while !region_member_parallelogram([0, 0], G₁ * 0.999999, G₂ * 0.999999, v₁)
                α = rand(-1:1)
                β = rand(-1:1)
                v₁ = v₀ + δB₂ + α * G₁ + β * G₂
            end
            position_1 = findfirst(x -> x == v₁, lattice_points[1:2*N_suc])
            if position_1 !== nothing
                push!(intralayer_connections, (x₁, position_1) => -Vppπ * exp(1im * [kx, ky] .* δB₂))
                push!(intralayer_connections, (position_1, x₁) => -Vppπ * exp(-1im * [kx, ky] .* δB₂))
            end

            # Third nearest neighbor
            v₁ = v₀ + δB₃
            while !region_member_parallelogram([0, 0], G₁ * 0.999999, G₂ * 0.999999, v₁)
                α = rand(-1:1)
                β = rand(-1:1)
                v₁ = v₀ + δB₃ + α * G₁ + β * G₂
            end
            position_1 = findfirst(x -> x == v₁, lattice_points[1:2*N_suc])
            if position_1 !== nothing
                push!(intralayer_connections, (x₁, position_1) => -Vppπ * exp(1im * [kx, ky] .* δB₃))
                push!(intralayer_connections, (position_1, x₁) => -Vppπ * exp(-1im * [kx, ky] .* δB₃))
            end

            # if verbose
            #     println("Processed x₁: $x₁")
            # end
        end
    elseif layer == 2
        # run over A sublattice of the twisted (2nd) layer
    end
end

function generate_interlayer_hopping_hamiltonian()
    return nothing
end
















# nrpointsuc = n^2 + n*m + m^2
# A1_points = Vector{Any}(undef, nrpointsuc)
# B1_points =  Vector{Any}(undef, nrpointsuc)
# bonds = Vector{Any}(undef, 3 * nrpointsuc)
# c₁ = 0
# c₂ = 0
# for (i, j) in product(-2*(n + m):2*(n + m), -2*(n + m):2*(n + m))
#     if region_member_parallelogram([0, 0], G₁ * 0.9999999, G₂ * 0.9999999, rA(i, j))
#         c₁ += 1
#         A1_points[c₁] = (i, j, rA(i, j))
#         c₂ += 1
#         bonds[c₂] = (rA(i, j), rA(i, j) + δB₁)
#         c₂ += 1
#         bonds[c₂] = (rA(i, j), rA(i, j) + δB₂)
#         c₂ += 1
#         bonds[c₂] = (rA(i, j), rA(i, j) + δB₃)
#     end
# end

# # resize A_points
# A1_points = [1:c₁]

# # reset counter
# c₁ = 0

# for (i, j) in product(-2*(n + m):2*(n + m), -2*(n + m):2*(n + m))
#     if region_member_parallelogram([0, 0], G₁, G₂, rB(i, j))
#         c₁ += 1
#         B1_points[c₁] = (i, j, rB(i, j))
#     end
# end

# # resize B_points 
# B1_points = B1_points[1:c₁]


# for (i, j) in product(-2*(n + m):2*(n + m), -2*(n + m):2*(n + m))
#     if region_member_parallelogram([0, 0], G₁ * 0.9999999, G₂ * 0.9999999, TrA(i, j))
#         c₁ += 1
#         A2_points[c₁] = (i, j, TrA(i, j))
#         c₂ += 1
#         Tbonds[c₂] = (TrA(i, j), TrA(i, j) + Tδ₁)
#         c₂ += 1
#         Tbonds[c₂] = (TrA(i, j), TrA(i, j) + Tδ₂)
#         c₂ += 1
#         Tbonds[c₂] = (TrA(i, j), TrA(i, j) + Tδ₃)
#     end
# end

# # resize A_points
# A2_points = A2_points[1:c₁]

# # reset counter
# c₁ = 0

# for (i, j) in product(-2*(n + m):2*(n + m), -2*(n + m):2*(n + m))
#     if region_member_parallelogram([0, 0], G₁ * 0.999999, G₂ * 0.999999, TrB(i, j))
#         c₁ += 1
#         B2_points[c₁] = (i, j, TrB(i, j))
#     end
# end

# # Resize B_points
# B2_points = B2_points[1:c₁]