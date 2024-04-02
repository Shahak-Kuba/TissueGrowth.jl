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
function sim1D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type)
    M = Int64(m*N)
    results = Vector{SimResults_t}(undef, 0)
    savetimes = LinRange(0, Tmax, 20)
    # creating 

    for ii in eachindex(D)
        @views kₛ = D[ii]*(η)/((l₀)^2)
        prob, p = SetupODEproblem1D(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir, dist_type)
        @time sol = solve(prob, RK4(), save_everystep = false, saveat=savetimes, dt=δt)
        #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
        push!(results, postSimulation1D(btype, sol, p))
        printInfo(ii,length(D),btype,N,kₛ,η,kf,M,D[ii])
    end

    return results

end 