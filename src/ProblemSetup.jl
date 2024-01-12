using DifferentialEquations

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
    u0 = ElasticArray{Float64}(undef,2,N)
    
    if btype == "circle"
        R = R₀ # to produce identical areas
        θ = collect(NodeDistribution(0.0,2*π,N+1,dist_type)) 
        pop!(θ)
        @views u0 .= [X(R,θ)'; Y(R,θ)'];
    elseif btype == "triangle"
        #R = √((2*π*R₀^2)/sin(π/3))
        R = √((π*R₀^2)/(√(3)*cos(π/6)^2))
        # calc verticies
        vertices = polygon_vertices(3, R, -π/2)
        # calc number of nodes per segment 
        w = Int64(N/3) + 1
        @views u0 .= position_vectors_polygon(vertices, w, dist_type)
    elseif btype == "square"
        #R = √(π*(R₀^2)) # to produce identical areas
        R = (R₀√(2π))/2
        # calc verticies
        vertices = polygon_vertices(4, R, -π/4)
        # calc number of nodes per segment 
        w = Int64(N/4) + 1
        @views u0 .= position_vectors_polygon(vertices, w, dist_type)
    elseif btype == "hex"
        R = √((2/3√3)*π*(R₀^2)) # to produce identical areas
        # calc verticies
        vertices = polygon_vertices(6, R, 0)
        # calc number of nodes per segment 
        w = Int64(N/6) + 1
        @views u0 .= position_vectors_polygon(vertices, w, dist_type)
    elseif btype == "star"
        star_points = 5
        Rotation_Angle = pi/2
        rotation_angle = Rotation_Angle + pi/star_points
        vertices = StarVerticies(star_points, R₀, Rotation_Angle, rotation_angle)
        w = Int64(N/(2star_points)) + 1
        u0 .= position_vectors_polygon(vertices, w, dist_type)
    elseif btype == "cross"
        side_length = √((π*R₀^2)/5)
        offset = side_length/2
        vertices = CrossVertecies(side_length, offset)
        w = Int64(N/12) + 1
        @views u0 .= position_vectors_polygon(vertices, w, dist_type)
    end
    return u0
end

"""
    SetupODEproblem1D(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir, dist_type)

Create and configure a 1D ODE problem for mechanical relaxation simulations.

This function sets up the initial conditions and parameters for a 1D ODE problem based on the specified boundary type and other physical parameters.

# Arguments
- `btype`: Boundary type (e.g., 'circle', 'triangle').
- `M`, `m`: Number of springs/particles in the system.
- `R₀`: Initial radius for circular boundary problems.
- `kₛ`, `η`: Prescaled mechanical relaxation coefficients.
- `kf`: Tissue production rate per cell per unit time.
- `l₀`: Resting spring length.
- `δt`: Euler timestep size.
- `Tmax`: End of simulation time.
- `growth_dir`: Direction of tissue growth ('inward' or 'outward').
- `dist_type`: Distribution type for node placement.

# Returns
Configured ODEProblem instance and parameters for the 1D simulation.
"""
function SetupODEproblem1D(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,growth_dir,dist_type)
    l₀ = l₀/m
    kₛ = kₛ*m
    kf = kf/m
    η = η/m
    # setting up initial conditions
    u0 = u0SetUp(btype,R₀,M,dist_type)
    # solving ODE problem
    p = (M,kₛ,η,kf,l₀,δt,growth_dir)
    tspan = (0.0,Tmax)
    return ODEProblem(ODE_fnc_1D!,u0,tspan,p), p
end


"""
    SetupODEproblem2D(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir, prolif, death, embed, α, β, γ, dist_type)

Set up and configure a 2D ODE problem for mechanical relaxation simulations in tissue growth.

This function initializes the conditions and parameters for a 2D ODE problem based on the specified boundary type, physical parameters, and cell behaviors. It then constructs an ODEProblem object, ready for solving with DifferentialEquations.jl.

# Arguments
- `btype`: Boundary type (e.g., 'circle', 'triangle').
- `M`: Total number of springs along the interface.
- `m`: Number of springs per cell.
- `R₀`: Initial radius or characteristic length of the shape.
- `kₛ`: Spring stiffness coefficient.
- `η`: Viscosity or damping coefficient.
- `kf`: Tissue production rate per cell.
- `l₀`: Resting length of the spring per cell.
- `δt`: Time step for the numerical integration.
- `Tmax`: Total simulation time.
- `growth_dir`: Direction of tissue growth ('inward' or 'outward').
- `prolif`: Boolean flag for cell proliferation.
- `death`: Boolean flag for cell death.
- `embed`: Boolean flag for cell embedding.
- `α, β, γ`: Parameters for the cell behaviors.
- `dist_type`: Distribution type for node placement.

# Returns
- `ODEProblem`: An ODE problem instance set up with the specified parameters and initial conditions.
- `p`: A tuple containing the parameters used in setting up the ODE problem.
"""
function SetupODEproblem2D(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,growth_dir,prolif,death,embed,α,β,γ,dist_type)
    l₀ = l₀/m
    kₛ = kₛ*m
    η = η/m
    kf = kf/m
    u0 = u0SetUp(btype,R₀,M,dist_type)
    #plotInitialCondition(u0)
    # solving ODE problem
    p = (m,kₛ,η,kf,l₀,δt,growth_dir,prolif,death,embed,α,β,γ)
    tspan = (0.0,Tmax)
    return ODEProblem(ODE_fnc_2D!,u0,tspan,p), p
end