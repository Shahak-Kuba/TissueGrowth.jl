using TissueGrowth
#include("Run/Comparison_Sims/ComparisonSimulation.jl")
include("ComparisonSimulation.jl")

# Shared variables
R₀ = 1.05#1.2694265629824517;
D_array = [0.0001, 0.01, 1];
kf = 0.006;
growth_dir = "inward";
Tmax = 26.0; # days
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
xbound = 1.1
ybound = 1.1
Cbar_min = 2
Cbar_max = 5
idx = 3
f2 = DiscVSContShape_plot(Discrete_Solution[idx], m, Continuum_Solution[idx], xbound, ybound, cmap, Cbar_min, Cbar_max)
save("Disc_VS_Cont_D_High.png", f2)