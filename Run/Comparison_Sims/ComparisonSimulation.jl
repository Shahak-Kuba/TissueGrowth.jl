include("ComparisonPlottingFncs.jl")

function ComparisonSim(N,m1,m2,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type, prolif, death, embed, α, β, γ, event_δt, seed, A)
    ### Discrete Simulation
    
    # simulation with m1 cells
    sols2D_m1, z, c = TissueGrowth.sim2D(N,m1,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,[btype],dist_type,
    prolif, death, embed, α, β, γ, event_δt, seed, 11);

    # simulation with m2 cells
    sols2D_m2, z, c = TissueGrowth.sim2D(N,m2,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,[btype],dist_type,
    prolif, death, embed, α, β, γ, event_δt, seed, 11);


    ### Continuum Simulation
    ρ₀ = sols2D_m1[1][1].Density[1][1];

    #using FVM for low diffusivity and FD for mid-high diffusivity
    if D >= 0.005
    θ_cont,R_cont,ρ_cont = TissueGrowth.FD_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,R₀,btype,growth_dir);
    else 
    θ_cont,R_cont,ρ_cont = TissueGrowth.FVM_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,R₀,btype, growth_dir);
    end


    # plotting
    Discrete_Solution_m1 = sols2D_m1[1][1];
    Discrete_Solution_m2 = sols2D_m2[1][1];
    Continuum_Solution = (θ_cont,R_cont,ρ_cont);
    indicies = [1,3,5,7,9,11]
    num_cols = 2
    f_results = TissueGrowth.DiscVSContDensity_plot_all(Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies, num_cols)    
    return Discrete_Solution_m1, Discrete_Solution_m2, Continuum_Solution, f_results
end