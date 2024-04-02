

"""
    X(R,θ), Y(R,θ), Xₜ(R,T), Yₜ(R,T), ...

Define boundary functions for various shapes like circles, triangles, squares, and hexagons.

These functions calculate the x and y coordinates for points on the boundary of these shapes based on given parameters.

# Arguments
- `R`: Radius or a characteristic length scale of the shape.
- `θ` or `T`: Parametric variable(s) used to define the position on the boundary.

# Returns
The x and y coordinates of a point on the boundary of the shape.
"""
# Circular Boundary
X(R,θ) = R*cos(θ);
Y(R,θ) = R*sin(θ);

# Triangluar Boundary
Xₜ(R,T) = ((0<=T) & (T<=1)) * (R/2-T*R) + 
    ((1<T) & (T<=2)) * (-R/2 + (T-1)*R/2) + 
    ((2<T) & (T<=3)) * ((T-2)*R/2) 

Yₜ(R,T) = ((0<=T) & (T<=1)) * ((R*sin(π/3))/2) + 
    ((1<T) & (T<=2)) * ((R*sin(π/3))/2 - (T-1)*(R*sin(π/3))) + 
    ((2<T) & (T<=3)) * (-(R*sin(π/3))/2 + (T-2)*(R*sin(π/3)))

# Square Boundary
Xₛ(R,T) = ((0<=T) & (T<=1)) * (R*T-R/2) + 
    ((1<T) & (T<=2)) * (R/2) + 
    ((2<T) & (T<=3)) * (R/2 - R*(T-2)) + 
    ((3<T) & (T<=4)) * (-R/2);
Yₛ(R,T) = ((0<=T) & (T<=1)) * (-R/2) + 
    ((1<T) & (T<=2)) * (-R/2 + R*(T-1)) + 
    ((2<T) & (T<=3)) * (R/2) + 
    ((3<T) & (T<=4)) * (R/2 - R*(T-3));;

# Hexagon Boundary
Vertex(R) = [R R*cos(pi/3) R*cos(2*pi/3) R*cos(pi) R*cos(4*pi/3) R*cos(5*pi/3) R*cos(2*pi); 0 R*sin(pi/3) R*sin(2*pi/3) R*sin(pi) R*sin(4*pi/3) R*sin(5*pi/3) R*sin(2*pi)];

Xₕ(R,T) = ((0<=T) & (T<=1)) * (Vertex(R)[1,1] + 
    T *(Vertex(R)[1,2] - Vertex(R)[1,1])) + 
    ((1<T) & (T<=2)) * (Vertex(R)[1,2] + 
    (T-1)*(Vertex(R)[1,3] - Vertex(R)[1,2]))+ 
    ((2<T) & (T<=3)) * (Vertex(R)[1,3] + 
    (T-2)*(Vertex(R)[1,4] - Vertex(R)[1,3]))+ 
    ((3<T) & (T<=4)) * (Vertex(R)[1,4] + 
    (T-3)*(Vertex(R)[1,5] - Vertex(R)[1,4]))+ 
    ((4<T) & (T<=5)) * (Vertex(R)[1,5] + 
    (T-4)*(Vertex(R)[1,6] - Vertex(R)[1,5]))+ 
    ((5<T) & (T<6))  * (Vertex(R)[1,6] + 
    (T-5)*(Vertex(R)[1,7] - Vertex(R)[1,6]))

Yₕ(R,T) = ((0<=T) & (T<=1)) * (Vertex(R)[2,1] + 
    T*(Vertex(R)[2,2] - Vertex(R)[2,1]))+ 
    ((1<T) & (T<=2)) * (Vertex(R)[2,2] + 
    (T-1)*(Vertex(R)[2,3] - Vertex(R)[2,2]))+ 
    ((2<T) & (T<=3)) * (Vertex(R)[2,3] + 
    (T-2)*(Vertex(R)[2,4] - Vertex(R)[2,3]))+ 
    ((3<T) & (T<=4)) * (Vertex(R)[2,4] + 
    (T-3)*(Vertex(R)[2,5] - Vertex(R)[2,4]))+ 
    ((4<T) & (T<=5)) * (Vertex(R)[2,5] + 
    (T-4)*(Vertex(R)[2,6] - Vertex(R)[2,5]))+ 
    ((5<T) & (T<6))  * (Vertex(R)[2,6] + 
    (T-5)*(Vertex(R)[2,7] - Vertex(R)[2,6]))

Xᵩ(T) = T
Yᵩ(T) = 2 + 0.5*cos(3*T)

##########################################################################################################################################################################


"""
    initial_pos_1D(u0, N, η, kf, l₀), initial_pos_2D(u0, N, η, kf, l₀)

Calculate the initial positions for 1D simulations.

These functions solve an ODE problem to determine the initial positions of particles or cells for the given simulation parameters.

# Arguments
- `u0`: Initial state vector.
- `N`: Number of particles or cells.
- `η`: Viscosity or a related parameter.
- `kf`: Growth factor or a related parameter.
- `l₀`: Characteristic length scale.

# Returns
The calculated initial positions for the simulation.
"""
function initial_pos_1D(u0,N,η,kf,l₀)
    kₛ = 2*N
    η = η/N
    kf = kf/N
    Tmax = 20;
    δt = 0.01
    p = (N,kₛ,η,kf,l₀,δt)
    tspan = (0.0,Tmax)
    prob = ODEProblem(ODE_fnc_1D_init!,u0,tspan,p)
    savetimes = LinRange(0, Tmax, 2)
    init_pos = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
    return init_pos.u[2]
end

"""
    initial_pos_2D(u0, N, η, kf, l₀)

Calculate the initial positions for 2D simulations.

These functions solve an ODE problem to determine the initial positions of particles or cells for the given simulation parameters.

# Arguments
- `u0`: Initial state vector.
- `N`: Number of particles or cells.
- `η`: Viscosity or a related parameter.
- `kf`: Growth factor or a related parameter.
- `l₀`: Characteristic length scale.

# Returns
The calculated initial positions for the simulation.
"""
function initial_pos_2D(u0,N,η,kf,l₀)
    kₛ = 1
    η = η/N
    kf = kf/N
    Tmax = 10;
    δt = 0.0005
    growth_dir = ""
    p = (N,kₛ,η,kf,l₀,δt,growth_dir)
    tspan = (0.0,Tmax)
    prob = ODEProblem(ODE_fnc_2D_init!,u0,tspan,p)
    savetimes = LinRange(0, Tmax, 2)
    init_pos = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
    return init_pos.u[2]
end

"""
    u0SetUp(btype, R₀, N, dist_type)

Set up initial conditions for simulations based on the boundary type and distribution.

This function initializes the positions of particles or cells based on the specified boundary type and distribution.

# Arguments
- `btype`: Type of boundary (e.g., 'circle', 'triangle').
- `R₀`: Initial radius or characteristic length.
- `N`: Number of points or particles.
- `dist_type`: Type of distribution for the points.

# Returns
An array of initial positions.
"""
function u0SetUp(btype,R₀,N,dist_type)
    # setting up initial conditions
    #θ = collect(LinRange(0.0, 2*π, N+1))  # just use collect(θ) to convert into a vector
    θ = collect(NodeDistribution(0.0,2*π,N+1,dist_type)) 
    pop!(θ)
    #u0 = Array{Float64}(undef,2,N)
    u0 = ElasticArray{Float64}(undef,2,N)
    for i in 1:N
        if btype == "circle"
            R = R₀ # to produce identical areas
            @views u0[:,i] .= [X(R,θ[i]), Y(R,θ[i])];
        elseif btype == "triangle"
            R = √((2*π*R₀^2)/sin(π/3))
            @views u0[:,i] .= [Xₜ(R,θ[i]*3/(2*π)), Yₜ(R,θ[i]*3/(2*π))]
        elseif btype == "square"
            R = √(π*(R₀^2)) # to produce identical areas
            @views u0[:,i] .= [Xₛ(R,θ[i]*2/pi), Yₛ(R,θ[i]*2/pi)]
        elseif btype == "hex"
            R = √((2/3√3)*π*(R₀^2)) # to produce identical areas
            @views u0[:,i] .= [Xₕ(R,θ[i]*3/pi), Yₕ(R,θ[i]*3/pi)]
        elseif btype == "SineWave" # length from (0→2π) ≈ 8.984
            θ = collect(LinRange(0.0, 2*π, N+1))  # just use collect(θ) to convert into a vector
            #pop!(θ)
            @views u0[:,i] .= [Xᵩ(θ[i]), Yᵩ(θ[i])];
        elseif btype == "star"
            star_points = 5
            r₀ = R₀/2 
            Rotation_Angle = pi/2
            rotation_angle = Rotation_Angle + pi/star_points
            verts = StarVerticies(star_points, R₀, Rotation_Angle, r₀, rotation_angle)
            u0 = interpolate_vertices(verts, Int( round(N/(2*star_points))))'
        elseif btype == "cross"
            side_length = √((π*R₀^2)/5)
            offset = side_length/2
            verts = CrossVertecies(side_length, offset)
            u0 = interpolate_vertices(verts, Int( round(N/12)))'
        end
    end

    return u0
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

"""
    nonLinearRange(start, stop, length, dist_type)

Generate a non-linear range of values between `start` and `stop`.

This function creates a range of values with various non-linear distributions like exponential, sinusoidal, etc.

# Arguments
- `start`: Starting value of the range.
- `stop`: Ending value of the range.
- `length`: Number of points in the range.
- `dist_type`: Type of non-linear distribution.

# Returns
A range of non-linearly distributed values.
"""
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


##############################################################################################

# linear interpolation for custom shapes

function interpolate_vertices(vertices, N)
    # Ensure at least two vertices are given
    if length(vertices) < 2
        error("At least two vertices are needed.")
    end

    interpolated_x = Float64[]
    interpolated_y = Float64[]

    for i in 1:size(vertices,1)-1
        x1, y1 = vertices[i,:]
        x2, y2 = vertices[i+1,:]

        # Calculate step sizes for x and y coordinates
        step_x = (x2 - x1) / N
        step_y = (y2 - y1) / N

        for j in 0:N-1
            x = x1 + j * step_x
            y = y1 + j * step_y
            push!(interpolated_x, x)
            push!(interpolated_y, y)
        end
    end

    return hcat([interpolated_x, interpolated_y]...)
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

function StarVerticies(N, R, Rotation_Angle, r, rotation_angle)
    # finding R such that areas will match with formula A = 2N(0.5*R₀*rₒ*sin(θ))
    R₀ = √((4*π*(R^2))/(2*N*sin(π/N)))
    rₒ = R₀/2
    # empty vector
    StarVerts = Vector{Float64}[]
    # generating verticies for outside polygon
    VERTS = regular_polygon_vertices(N, R₀, Rotation_Angle)
    # generating verticies for inside polygon
    verts = regular_polygon_vertices(N, rₒ, rotation_angle)

    # combining the two
    for i in 1:N
        push!(StarVerts, VERTS[i,:])
        push!(StarVerts, verts[i,:])
    end
    push!(StarVerts, VERTS[1,:])

    return hcat(StarVerts...)'
end


# Cross Interface

function CrossVertecies(side_length, offset)
    CrossVerts = Vector{Float64}[]
    push!(CrossVerts,[side_length + offset, offset])
    push!(CrossVerts,[offset, offset])
    push!(CrossVerts,[offset, offset + side_length])
    push!(CrossVerts,[-offset, offset + side_length])
    push!(CrossVerts,[-offset, offset])
    push!(CrossVerts,[-offset - side_length, offset])
    push!(CrossVerts,[-offset - side_length, -offset])
    push!(CrossVerts,[-offset, -offset])
    push!(CrossVerts,[-offset, -offset - side_length])
    push!(CrossVerts,[offset, -offset - side_length])
    push!(CrossVerts,[offset, -offset])
    push!(CrossVerts,[offset + side_length, -offset])
    push!(CrossVerts,[side_length + offset, offset])

    return hcat(CrossVerts...)'
end

