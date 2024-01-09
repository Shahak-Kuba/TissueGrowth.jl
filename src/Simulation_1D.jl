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