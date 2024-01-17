"""
    sim1D()

Execute a series of 1D mechanical relaxation simulations.

This function sets up and runs a series of 1D simulations for different stiffness values. It iterates over an array of stiffness coefficients, sets up the corresponding ODE problem for each case, solves it, and collects the results.

# Simulation Parameters
- `N`: Number of cells in the simulation.
- `m`: Number of springs per cell.
- `M`: Total number of springs along the interface.
- `R₀`: Radius or characteristic length of the initial shape.
- `D`: Array of diffuision coefficients used to calculate cell stiffness.
- `l₀`: Resting length of the spring per cell.
- `kf`: Tissue production rate per cell.
- `η`: Viscosity or damping coefficient per cell.
- `Tmax`: Total simulation time (in days).
- `δt`: Time step for the Euler integration.
- `growth_dir`: Growth direction for the simulation ('1D').
- `btype`: Type of boundary condition ('SineWave').
- `savetimes`: Time points at which to save the simulation results.

# Returns
A vector of `SimResults_t` objects, each representing the simulation results for a different stiffness coefficient.
"""
function sim1D()
    # setting up simulation parameters
    N = 77 # number of cells
    m = 2 # number of springs per cell
    M = m*N # total number of springs along the interface
    R₀ = 1  # shape radius
    #kₛ_Array = [0.1, 1, 5, 10]
    D = [0.001, 0.075, 0.15, 1] # of cell 
    l₀ = 1 # of cell 
    kf = 0.01*m # of cell = kf¹
    η = 1 # of cell = η¹
    Tmax = 25 # days
    δt = 0.0001
    growth_dir = "1D"
    #btypes = ["circle", "triangle", "square", "hex"]
    btype = "SineWave"
    savetimes = LinRange(0, Tmax, 30)

    #sol_array = Array{ODESolution}(undef,length(btypes));
    results = Vector{SimResults_t}(undef, 0)
    # creating 

    for ii in eachindex(D)
        @views kₛ = D[ii]*(η)/((l₀)^2)
        prob, p = SetupODEproblem1D(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir)
        @time sol = solve(prob, Euler(), save_everystep = false, saveat=savetimes, dt=δt)
        #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
        push!(results, postSimulation1D(btype, sol, p))
        printInfo(ii,length(D),btype,N,kₛ,η,kf,M,D[ii])
    end

    return results

end 