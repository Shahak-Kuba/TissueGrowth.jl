using BenchmarkTools

"""
    sim2D()

Execute a series of 2D mechanical relaxation simulations.

This function sets up and runs 2D simulations for different stiffness coefficients and boundary types. It iterates over arrays of stiffness coefficients and boundary types, sets up the corresponding ODE problem for each case, solves it using a specified numerical method with periodic callbacks, and collects the results.

# Simulation Parameters
- `N`: Number of cells in the simulation.
- `m`: Number of springs per cell.
- `M`: Total number of springs along the interface.
- `R₀`: Radius or characteristic length of the initial shape.
- `D`: Array of diffuision coefficients used to calculate cell stiffness.
- `l₀`: Resting length of the spring per cell.
- `kf`: Tissue production rate per cell.
- `η`: Viscosity or damping coefficient per cell.
- `growth_dir`: Direction of tissue growth ('inward' or 'outward').
- `Tmax`: Total simulation time (in days).
- `δt`: Time step for the numerical integration.
- `btypes`: Types of boundary conditions (e.g., 'circle', 'triangle').
- `dist_type`: Distribution type for node placement (e.g., 'Linear', 'sigmoid').
- `prolif`, `death`, `embed`: Boolean flags indicating cell behaviors.
- `α`, `β`, `γ`: Parameters for cell behaviors.
- `event_δt`: Time interval for periodic callback events.
- `savetimes`: Time points at which to save the simulation results.

# Returns
A vector of vectors of `SimResults_t` objects. Each inner vector represents the simulation results for different boundary types under a specific stiffness coefficient.

# Example
```julia
all_results = sim2D()
```
"""
function sim2D()
    # setting up simulation parameters
    N = 180 # number of cells
    m = 1 # number of springs per cell
    M = Int(m*N) # total number of springs along the interface
    R₀ = 1  # shape radius
    D = [0.0001, 0.0075, 0.015]
    l₀ = 1
    kf = 0.0008
    η = 1
    growth_dir = "inward" # Options: "inward", "outward"
    Tmax = 21 # days
    δt = 0.01
    btypes = ["circle"]#["circle", "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
    dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

    ## Cell Behaviours
    prolif = true; death = true; embed = false;
       α = 0.01;        β = 0.1;      γ = 0.01;
    event_δt = 0.3

    savetimes = LinRange(0, Tmax, 8)

    all_results = Vector{Vector{SimResults_t}}(undef, 0)

    for jj in eachindex(D)
        @views kₛ = D[jj]*(η)/((l₀)^2)
        #sol_array = Array{ODESolution}(undef,length(btypes));
        results = Vector{SimResults_t}(undef, 0)
        # creating 

        for ii in eachindex(btypes)
            @views btype = btypes[ii]
            prob, p = SetupODEproblem2D(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,
                                        growth_dir,prolif,death,embed,α,β,γ,dist_type)
            cb = PeriodicCallback(affect!,event_δt; save_positions=(false, false))
            @time sol = solve(prob, RK4(), save_everystep = false, saveat=savetimes, dt=δt, dtmax = δt, callback = cb)
            push!(results, postSimulation2D(btype, sol, p))
            printInfo(ii,length(btypes),btype,N,kₛ*m,η/m,kf/m,M,D[jj])
        end
        push!(all_results,results)
    end

    return all_results

end
