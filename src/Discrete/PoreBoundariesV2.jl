using LinearAlgebra
using QuadGK
using Roots

# Circular Boundary
X(R,θ) = R.*cos.(θ);
Y(R,θ) = R.*sin.(θ);

# for all other regular regular polygons

function polygon_vertices(N, R, theta)
    # N: Number of sides
    # R: Radius
    # theta: Initial rotation angle in radians

    vertices = []
    for i in 0:(N-1)
        angle = theta + (2 * π * i / N)
        x = R * cos(angle)
        y = R * sin(angle)
        push!(vertices, (x, y))
    end
    return vertices
end

function interpolate_segment(p1, p2, w, dist_type)
    [((1 - t) .* p1 .+ t .* p2) for t in NodeDistribution(0, 1, w, dist_type)]
end

# for irregular polygons
function CrossVertecies(side_length, offset)
    CrossVerts = []
    push!(CrossVerts,(side_length + offset, offset))
    push!(CrossVerts,(offset, offset))
    push!(CrossVerts,(offset, offset + side_length))
    push!(CrossVerts,(-offset, offset + side_length))
    push!(CrossVerts,(-offset, offset))
    push!(CrossVerts,(-offset - side_length, offset))
    push!(CrossVerts,(-offset - side_length, -offset))
    push!(CrossVerts,(-offset, -offset))
    push!(CrossVerts,(-offset, -offset - side_length))
    push!(CrossVerts,(offset, -offset - side_length))
    push!(CrossVerts,(offset, -offset))
    push!(CrossVerts,(offset + side_length, -offset))
    #push!(CrossVerts,(side_length + offset, offset))

    #return hcat(CrossVerts...)'
    return CrossVerts
end

function StarVerticies(N, R, Rotation_Angle, rotation_angle)
    # finding R such that areas will match with formula A = 2N(0.5*R₀*rₒ*sin(θ))
    R₀ = √((4*π*(R^2))/(2*N*sin(π/N)))
    rₒ = R₀/2
    # empty vector
    StarVerts = []
    # generating verticies for outside polygon
    VERTS = regular_polygon_vertices(N, R₀, Rotation_Angle)
    # generating verticies for inside polygon
    verts = regular_polygon_vertices(N, rₒ, rotation_angle)

    # combining the two
    for i in 1:N
        push!(StarVerts, (VERTS[i,1], VERTS[i,2]))
        push!(StarVerts, (verts[i,1], verts[i,2]))
    end
    #push!(StarVerts, VERTS[1,:])

    #return hcat(StarVerts...)'
    return StarVerts
end

function regular_polygon_vertices(N, R, rotation_angle)
    vertices = Vector{Float64}[]

    for i in 0:N-1
        angle = 2π * i / N + rotation_angle
        x = R * cos(angle)
        y = R * sin(angle)
        push!(vertices, [x, y])
    end

    return hcat(vertices...)'
end

function position_vectors_polygon(vertices, N, dist_type)
    V = length(vertices)
    all_x = Float64[]
    all_y = Float64[]
    for i in 1:V
        # Handle wrap-around at the last vertex
        p1 = vertices[i]
        p2 = vertices[i % V + 1]
        segment_points = interpolate_segment(p1, p2, N, dist_type)

        for i in 1:(length(segment_points) - 1)
            point = segment_points[i]
            push!(all_x, point[1])
            push!(all_y, point[2])
        end
    end
    return [all_x'; all_y']
end

"""
    NodeDistribution(start, stop, length, type)

Generate a distribution of nodes between `start` and `stop` with a specified distribution type.

# Arguments
- `start`: Starting value of the range.
- `stop`: Ending value of the range.
- `length`: Number of points in the distribution.
- `type`: Type of distribution (e.g., 'Linear', 'exp', 'sine').

# Returns
A range or array of distributed values.
"""
function NodeDistribution(start,stop,length,type)
    if type == "Linear"
        return LinRange(start, stop, length)
    else
        return nonLinearRange(start, stop, length, type)
    end
end

function nonLinearRange(start, stop, length, dist_type)
    linear_range = LinRange(0, 1, length)

    # Applying different distribution types
    if dist_type == "exp"
        # Exponential scaling
        return start .+ (exp.(linear_range .* log(1 + stop - start)) .- 1)
    elseif dist_type == "sine"
        # Sinusoidal scaling
        return start .+ (sin.((π/2) .* linear_range) .* (stop - start))
    elseif dist_type == "cosine"
        # Cosine scaling
        return start .+ ((1 .- cos.((π/2) .* linear_range)) .* (stop - start))
    elseif dist_type == "quad"
        # Quadratic scaling
        return start .+ (linear_range .^ 2 .* (stop - start))
    elseif dist_type == "cubic"
        # Cubic scaling
        return start .+ (linear_range .^ 3 .* (stop - start))
    elseif dist_type == "sigmoid"
        # Sigmoid scaling
        linear_range = LinRange(-1, 1, length)
        k = 5; # slope steepness
        sigmoid_range = 1 ./ (1 .+ exp.(-k.*(linear_range)))
        return start .+ (sigmoid_range .* (stop - start))
    elseif dist_type == "2sigmoid"
        # Piecewise sigmoid scaling
        k1 = 10;  k2 = 10;
        x01 = 0.5;  x02 = 0.5;
        piecewise_sigmoid = [x < 0.5 ? 0.5 * (1 / (1 + exp(-k1 * (2x - x01)))) : 0.5 + 0.5 * (1 / (1 + exp(-k2 * (2x - 1 - x02)))) for x in linear_range]
        return start .+ (piecewise_sigmoid * (stop - start))
    elseif dist_type == "4sigmoid"
        # Parameters for the sigmoid functions
        k = 20
        # Adjust midpoints for the full sigmoid in the first and last segments
        x0 = [0.125, 0.375, 0.625, 0.875]

        # Piecewise sigmoid scaling with 4 segments
        piecewise_sigmoid = [if x < 0.25
                                (1 / (1 + exp(-k * (4x - x0[1])))) * 0.25
                             elseif x < 0.5
                                0.25 + (1 / (1 + exp(-k * (4x - 1 - x0[2])))) * 0.25
                             elseif x < 0.75
                                0.5 + (1 / (1 + exp(-k * (4x - 2 - x0[3])))) * 0.25
                             else
                                0.75 + (1 / (1 + exp(-k * (4x - 3 - x0[4])))) * 0.25
                             end for x in linear_range]

        return start .+ (piecewise_sigmoid .* (stop - start))
    else
        error("Unsupported distribution type")
    end
end

