using TissueGrowth
using Makie
include("../C-Sim_Code/Pore_Simulations/Discrete_Continuum_Comparison_Sims/ComparisonSimulation.jl")
include("../C-Sim_Code/Pore_Simulations/Discrete_Continuum_Comparison_Sims/ComparisonPlottingFncs.jl")

# Shared variables
R₀ = 56.41895835477563
D_array = [1, 100, 10000];
kf = 87.84;
growth_dir = "inward";
Tmax = 5; # days
btype = "square"; #Options: ["circle", "triangle", "square", "hex", "star","cross"]

# Discrete Simulation Variables
# set random seed number for reproducability 
seed = 88;

# setting up simulation parameters
N = 20; # number of cells
m = 10; # number of springs per cell
l₀ = 10.0;
η = 1.0 ;
δt = 0.01;
dist_type = "Linear"; #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        βv = 0.001;      γv = 0.01;
event_δt = δt;

# Continuum simulation variavbles
Av = 0.0;

# Generating results
Discrete_Solution, Continuum_Solution = ComparisonSim_Density(N,m,R₀,D_array,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type, 
                                                                            prolif, death, embed, α, βv, γv, event_δt, seed, Av);

cmap = :cool
xbound = 60
ybound = 60
Cbar_min = 0.05
Cbar_max = 0.15
f1 = TissueGrowth.DiscVSContShape_plot(Discrete_Solution[1], m, Continuum_Solution[1], xbound, ybound, cmap, Cbar_min, Cbar_max)
#save("Disc_VS_Cont_D_Low.png", f1)
f2 = TissueGrowth.DiscVSContShape_plot(Discrete_Solution[2], m, Continuum_Solution[2], xbound, ybound, cmap, Cbar_min, Cbar_max)
#save("Disc_VS_Cont_D_Mid.png", f2)
f3 = TissueGrowth.DiscVSContShape_plot(Discrete_Solution[3], m, Continuum_Solution[3], xbound, ybound, cmap, Cbar_min, Cbar_max)
#save("Disc_VS_Cont_D_High.png", f3)

f4 = TissueGrowth.DiscVSContShape_plot_all(Discrete_Solution,m,Continuum_Solution, xbound, ybound, cmap, Cbar_min, Cbar_max)
save("Disc_VS_Cont_D.png", f4)