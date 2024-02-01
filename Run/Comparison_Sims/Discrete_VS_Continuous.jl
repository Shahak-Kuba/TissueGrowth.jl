using TissueGrowth
#include("Run/Comparison_Sims/ComparisonSimulation.jl")
include("ComparisonSimulation.jl")

# Shared variables
R₀ = 1.05;
D = 0.001;
kf = 0.005;
growth_dir = "inward";
Tmax = 26.0; # days
btype = "square"; #Options: ["circle", "triangle", "square", "hex", "star","cross"]

# Discrete Simulation Variables
# set random seed number for reproducability 
seed = 88;

# setting up simulation parameters
N = 24; # number of cells
m1 = 1; # number of springs per cell
m2 = 4;
l₀ = 1.0;
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
Discrete_Solution_m1, Discrete_Solution_m2, Continuum_Solution, f_results = ComparisonSim(N,m1,m2,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type, 
                                                                            prolif, death, embed, α, βv, γv, event_δt, seed, Av);

indicies = [1,3,5,7,9,11]
num_cols = 3
f1 = DiscVSContDensity_plot_all(Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies, num_cols)

cmap = :jet
xbound = 1.1
ybound = 1.1
Cbar_min = 0
Cbar_max = 10
f2 = DiscVSContShape_plot(Discrete_Solution_m2, m2, Continuum_Solution, xbound, ybound, cmap, Cbar_min, Cbar_max)
