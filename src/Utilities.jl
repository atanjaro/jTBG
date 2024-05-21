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
function rotate_vector(vector, θ)
    # Convert vector to polar coordinates
    magnitude = abs(vector)
    current_angle = atan(imag(vector), real(vector))

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


"""
    in_parallelogram(p1, p2, p)

Checks whether the point 'p' lies in the parallelogram formed by 'p₁' and 'p₂'
"""
function in_parallelogram(p₁, p₂, p)
    v₁ = p₂ .- p₁
    v₂ = p .- p₁
    
    # Calculate the cross product
    cross_product = la.cross(v₁, v₂)
    
    # Check if the cross product is approximately zero
    return la.norm(cross_product) < 1e-10
end

function parallelogram_corner(a, b)
    c = a + b
    d = c + a
    Polygon([a, b, c, d])
end


function region_member(G1, G2, rA, i, j)
    parallelogram = parallelogram_corner((0, 0), (G1 * 0.9999999, G2 * 0.9999999))
    rA_point = rA(i, j)
    rA_point in parallelogram
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