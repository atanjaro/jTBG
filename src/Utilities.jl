"""
    angle_between_vectors( v1:Vector{T}, v2::Vector{T} )

Finds the angle between vectors 'v₁' and 'v₂' in radians.
"""
function angle_between_vectors(v₁,v₂)
    dot_product = la.dot(v₁, v₂)
    norm_v₁ = la.norm(v₁)
    norm_v₂ = la.norm(v₂)
    
    # Avoid division by zero
    if norm_v₁ == 0 || norm_v₂ == 0
        error("One of the vectors has zero length.")
    end
    
    cos_angle = dot_product / (norm_v₁ * norm_v₂)
    angle_in_radians = acos(cos_angle)
    
    return angle_in_radians
end


"""
    rotate_vector( vector::Complex{T}, θ::Real )

Rotates a vector by an angle specified in radians.
"""
# function rotate_vector(vector, θ)
#     # Convert vector to polar coordinates
#     magnitude = la.norm(vector)
#     current_angle = atan(imag(vector), real(vector))

#     # Apply rotation to the angle component
#     new_angle = current_angle + θ

#     # Convert back to Cartesian coordinates
#     rotated_vector = magnitude * exp(im * new_angle)

#     return rotated_vector
# end

function rotate_vector(vector, θ::Float64)
    if isa(vector, Complex)
        # If vector is a complex number
        magnitude = norm(vector)
        current_angle = angle(vector)
    elseif isa(vector, AbstractVector) && length(vector) == 2
        # If vector is a 2D Cartesian coordinate
        real_part, imag_part = vector
        magnitude = la.norm(vector)
        current_angle = atan(imag_part, real_part)
    else
        error("Unsupported vector format. Provide a complex number or a 2D Cartesian coordinate.")
    end

    # Apply rotation to the angle component
    new_angle = current_angle + θ

    # Convert back to Cartesian coordinates
    rotated_vector = magnitude * exp(im * new_angle)

    return rotated_vector
end



"""
    rotation_matrix( θ::Real )

Returns a rotation matrix given a specified angle 'θ' (in radians).
"""
function rotation_matrix(θ)
    cosθ = cos(θ)
    sinθ = sin(θ)
    return [cosθ -sinθ; sinθ cosθ]
end


# Function to check if point is in the parallelogram
function region_member_parallelogram(origin, v1, v2, point)
    relative_point = point .- origin
    dot_v1_v1 = la.dot(v1, v1)
    dot_v1_v2 = la.dot(v1, v2)
    dot_v2_v2 = la.dot(v2, v2)
    dot_rp_v1 = la.dot(relative_point, v1)
    dot_rp_v2 = la.dot(relative_point, v2)
    det = dot_v1_v1 * dot_v2_v2 - dot_v1_v2 * dot_v1_v2
    alpha = (dot_v2_v2 * dot_rp_v1 - dot_v1_v2 * dot_rp_v2) / det
    beta = (dot_v1_v1 * dot_rp_v2 - dot_v1_v2 * dot_rp_v1) / det
    return 0 <= alpha <= 1 && 0 <= beta <= 1
end


"""
    get_kq_points( k_points::Vector{SVector{2, Float64}}, q_point::Vector{Float64} )

Given a set of k-points along cuts in the FBZ and a valid q-point, computes (k+q)-points, where q is the phonon momentum. 

"""
function get_kq_points(k_points, q_point)        
    kq_points = []                      
    for i in k_points
        push!(kq_points,i.+q_point)
    end

    return kq_points
end


"""
    get_all_kq_points( k_points::Vector{SVector{2, Float64}}, q_points::Vector{Float64} )

Given a set of k-points along cuts in the FBZ and a valid q-point, computes all (k+q)-points, where q is the phonon momentum. 

"""
function get_all_kq_points(k_points, q_points)
    all_kq_points = []
    for q_point in q_points
        points = get_kq_points(k_points, q_point)
        push!(all_kq_points, points)
    end

    return all_kq_points        # creates an array of arrays of q points, where each sub-array corresponds to one particular q-points
end



"""
    get_positions(a₁, a₂, Nk)

Obtains position vectors for a finite lattice given the number of k-points. 
"""
function get_positions(a₁, a₂, Nk)
    rs = []
    for i in 0:(Nk)-1
        for j in 0:Nk-1
            r = i* a₁ + j* a₂
            push!(rs,r)
        end
    end

    return rs
end


# various convenience functions
rA(i, j) = i*a₁ + j*a₂
rB(i, j) = (rA(i,j) + δB₁)
TrA(i, j) = i * Ta₁ + j * Ta₂ + TδA
TrB(i, j) = i * Ta₁ + j * Ta₂ + TδB