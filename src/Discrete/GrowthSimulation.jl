"""
    sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,prolif,death,embed,α,β,γ,event_δt,seed, NumSaveTimePoints)

Execute a series of 2D mechanical relaxation simulations.

This function sets up and runs 2D simulations for different stiffness coefficients and boundary types. It iterates over arrays of stiffness coefficients and boundary types, sets up the corresponding ODE problem for each case, solves it using a specified numerical method with periodic callbacks, and collects the results.

# Simulation Parameters
- `N`: Number of cells in the simulation.
- `m`: Number of springs per cell.
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
- `seed`: Seed number for reproducability with stochastic simulations.

# Returns
A vector of vectors of `SimResults_t` objects. Each inner vector represents the simulation results for different boundary types under a specific stiffness coefficient.

# Example
```julia
all_results = sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,
                        prolif, death, embed, α, β, γ, event_δt, seed, NumSaveTimePoints);
```
"""
function GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,restoring_force,dist_type, 
                prolif, death, embed, β, γ, Ot, event_δt, seed, NumSaveTimePoints)

    #Set_Random_Seed(seed)
    M = Int(m*N) # total number of springs along the interface
    savetimes = LinRange(0, Tmax, NumSaveTimePoints)
    
    # for cell embedment 
    global embedded_cells = []
    global cell_embedment_times = []
    embedded_cells_count = []

    #all_results = Vector{Vector{SimResults_t}}(undef, 0)

    event_cb = PeriodicCallback(event_affect!,event_δt; save_positions=(false, false))

    # saving callback 
    saved_values = SavedValues(Float64, Float64)
    save_cb = SavingCallback(store_embedded_cell_count, saved_values, saveat=savetimes)
    
    # generating a set of callbacks
    cbs = CallbackSet(event_cb,save_cb)

    if D > 0.0 # Simulations for constant diffusion (Nonlinear restoring force) where user specifies diffusivity
        kₛ = D*(η)/((l₀)^2)
        results = Vector{SimResults_t}(undef, 0)

        for ii in eachindex(btypes)
            @views btype = btypes[ii]
            prob, p = SetupODEproblem(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,
                                        growth_dir,domain_type,prolif,death,embed,β,γ,Ot,dist_type,restoring_force)
            Set_Random_Seed(seed)                            
            @time sol = solve(prob, RK4(), save_everystep = false, saveat=savetimes, dt=δt, dtmax = δt, callback = cbs)
            push!(results, postSimulation(btype, sol, p))
            push!(embedded_cells_count, floor.(saved_values.saveval))
            printInfo(ii,length(btypes),btype,N,kₛ*m,η/m,kf/m,M,D)
        end

    else # simulations where user specifies spring stiffness (used for a general restoring force)
        results = Vector{SimResults_t}(undef, 0)
        for ii in eachindex(btypes)
            @views btype = btypes[ii]
            prob, p = SetupODEproblem(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,
                                        growth_dir,domain_type,prolif,death,embed,β,γ,Ot,dist_type,restoring_force)
            Set_Random_Seed(seed)
            @time sol = solve(prob, RK4(), save_everystep = false, saveat=savetimes, dt=δt, dtmax = δt, callback = cbs)
            push!(results, postSimulation(btype, sol, p))
            push!(embedded_cells_count, floor.(saved_values.saveval))
            printInfo(ii,length(btypes),btype,N,kₛ*m,η/m,kf/m,M,D)
        end
    end

    return results, convert_matrix(hcat(embedded_cells...),m+1), embedded_cells_count, cell_embedment_times

end


# Growth simulation for when cell denisty limit is applied
function GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,restoring_force,dist_type, 
    prolif, death, embed, β, γ, Ot, event_δt, seed, NumSaveTimePoints, q_lim)

        #Set_Random_Seed(seed)
        M = Int(m*N) # total number of springs along the interface
        savetimes = LinRange(0, Tmax, NumSaveTimePoints)

        # for cell embedment 
        global embedded_cells = []
        global cell_embedment_times = []
        embedded_cells_count = []

        all_results = Vector{Vector{SimResults_t}}(undef, 0)

        event_cb = PeriodicCallback(event_affect!,event_δt; save_positions=(false, false))
        q_lim_cb = DiscreteCallback(density_lim_condition, density_lim_affect!)
        # saving callback 
        saved_values = SavedValues(Float64, Float64)
        save_cb = SavingCallback(store_embedded_cell_count, saved_values, saveat=savetimes)

        # generating a set of callbacks
        cbs = CallbackSet(event_cb,q_lim_cb,save_cb)

        if D > 0.0 # Simulations for constant diffusion (Nonlinear restoring force) where user specifies diffusivity
            kₛ = D*(η)/((l₀)^2)
            results = Vector{SimResults_t}(undef, 0)
    
            for ii in eachindex(btypes)
                @views btype = btypes[ii]
                prob, p = SetupODEproblem(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,
                                            growth_dir,domain_type,prolif,death,embed,β,γ,Ot,dist_type,restoring_force)
                Set_Random_Seed(seed)
                @time sol = solve(prob, RK4(), save_everystep = false, saveat=savetimes, dt=δt, dtmax = δt, callback = cbs)
                push!(results, postSimulation(btype, sol, p))
                push!(embedded_cells_count, floor.(saved_values.saveval))
                printInfo(ii,length(btypes),btype,N,kₛ*m,η/m,kf/m,M,D)
            end
    
        else # simulations where user specifies spring stiffness (used for a general restoring force)
            results = Vector{SimResults_t}(undef, 0)
            for ii in eachindex(btypes)
                @views btype = btypes[ii]
                prob, p = SetupODEproblem(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,
                                            growth_dir,domain_type,prolif,death,embed,β,γ,Ot,dist_type,restoring_force)
                Set_Random_Seed(seed)
                @time sol = solve(prob, RK4(), save_everystep = false, saveat=savetimes, dt=δt, dtmax = δt, callback = cbs)
                push!(results, postSimulation(btype, sol, p))
                push!(embedded_cells_count, floor.(saved_values.saveval))
                printInfo(ii,length(btypes),btype,N,kₛ*m,η/m,kf/m,M,D)
            end
        end

        return results, convert_matrix(hcat(embedded_cells...),m+1), embedded_cells_count, cell_embedment_times

end