using BenchmarkTools

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
