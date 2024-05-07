include("ComparisonPlottingFncs.jl")

function ComparisonSim(N,m1,m2,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type, prolif, death, embed, α, β, γ, event_δt, seed, Av)
    ### Discrete Simulation
    Λ = 10000
    
    # simulation with m1 cells
    sols2D_m1, z, c = TissueGrowth.GrowthSimulation(N,m1,R₀,D*Λ,0,l₀,kf,η,growth_dir,"2D",Tmax,δt,[btype],dist_type,
    prolif, death, embed, α, β, γ, event_δt, seed, 11);

    # simulation with m2 cells
    sols2D_m2, z, c = TissueGrowth.GrowthSimulation(N,m2,R₀,D*Λ,0,l₀,kf,η,growth_dir,"2D",Tmax,δt,[btype],dist_type,
    prolif, death, embed, α, β, γ, event_δt, seed, 11);


    ### Continuum Simulation
    ρ₀ = sols2D_m1[1].Density[1][1];

    #using FVM for low diffusivity and FD for mid-high diffusivity
    if D >= 0.005
        θ_cont,R_cont,ρ_cont = TissueGrowth.FD_SolveContinuumLim_Polar(D,kf,Av,ρ₀,Tmax,R₀,btype,growth_dir);
    else 
        θ_cont,R_cont,ρ_cont = TissueGrowth.FVM_SolveContinuumLim_Polar(D,kf,Av,ρ₀,Tmax,R₀,btype, growth_dir);
    end


    # plotting
    Discrete_Solution_m1 = sols2D_m1[1];
    Discrete_Solution_m2 = sols2D_m2[1];
    Continuum_Solution = (θ_cont,R_cont,ρ_cont);
    return Discrete_Solution_m1, Discrete_Solution_m2, Continuum_Solution
end

function ComparisonSim_Density(N,m,R₀,D_array,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type, prolif, death, embed, α, β, γ, event_δt, seed, Av)
    Discrete_Solution = [];
    Continuum_Solution = [];

    for D in D_array
   
        ### Discrete Simulation
        
        # simulation with m1 cells
        sol_Discrete, z, c = TissueGrowth.GrowthSimulation(N,m,R₀,D,0,l₀,kf,η,growth_dir,"2D",Tmax,δt,[btype],dist_type,
        prolif, death, embed, α, β, γ, event_δt, seed, 11);

        ### Continuum Simulation
        ρ₀ = sol_Discrete[1].Density[1][1];
        Λ = 100000
        D_cont = D / Λ

        #using FVM for low diffusivity and FD for mid-high diffusivity
        if D_cont >= 0.005
            θ_cont,R_cont,ρ_cont = TissueGrowth.FD_SolveContinuumLim_Polar(D,kf,Av,ρ₀,Tmax,R₀,btype,growth_dir);
        else 
            θ_cont,R_cont,ρ_cont = TissueGrowth.FVM_SolveContinuumLim_Polar(D,kf,Av,ρ₀,Tmax,R₀,btype, growth_dir);
        end


        # plotting
        push!(Discrete_Solution, sol_Discrete[1]);
        push!(Continuum_Solution,(θ_cont,R_cont,ρ_cont));
    end

    return Discrete_Solution, Continuum_Solution
end