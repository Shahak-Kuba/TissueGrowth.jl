function sim2D_δt()
    # setting up simulation parameters
    N = 180 # number of cells
    m = 1 # number of springs per cell
    M = Int(m*N) # total number of springs along the interface
    #N = 500
    R₀ = 1  # shape radius
    D = 0.0075
    l₀ = 1
    kf = 0.0008
    η = 1
    kₛ = D*(η)/((l₀)^2)
    Tmax = 20# days
    growth_dir = "inward" 
    δt_Array = [0.00001,0.0001,0.001,0.005]
     #["circle", "triangle", "square", "hex", "star","cross"]
    btypes = ["triangle"] #, "square", "hex", "star","cross"]
    dist_type = "Linear"
    savetimes = LinRange(0, Tmax, 30)
    rounded_savetimes = [round(x,digits=2) for x in savetimes]

    ## Cell Behaviours
    prolif = false; death = true; embed = false;
       α = 0.1;        β = 0.1;      γ = 0.01;
    event_δt = 0.1

    all_results = Vector{Vector{SimResults_t}}(undef, 0)

    for jj in eachindex(btypes)
        results = Vector{SimResults_t}(undef, 0)
        @views btype = btypes[jj]

        for ii in eachindex(δt_Array)
            @views δt = δt_Array[ii]
            # creating 
            prob, p = SetupODEproblem2D(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,
                                        growth_dir,prolif,death,embed,α,β,γ,dist_type)
            @time sol = solve(prob, Euler(), save_everystep = false, saveat=rounded_savetimes, dt=δt)
            #@btime sol = solve(prob, Euler(), saveat=savetimes, dt=δt)
            push!(results, postSimulation2D(btype, sol, p))
            printInfo(ii,length(btypes),btype,N,kₛ*m,η/m,kf/m,M,D)
        end

        push!(all_results,results)
    end

    return all_results

end