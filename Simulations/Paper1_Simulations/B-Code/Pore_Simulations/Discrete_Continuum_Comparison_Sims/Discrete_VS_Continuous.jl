using TissueGrowth
#include("Run/Comparison_Sims/ComparisonSimulation.jl")
include("ComparisonSimulation.jl")

# Shared variables
R₀ = 1.05#1.2694265629824517;
D = 0.0075;
kf = 0.006;
growth_dir = "inward";
Tmax = 22.0; # days
btype = "square"; #Options: ["circle", "triangle", "square", "hex", "star","cross"]

# Discrete Simulation Variables
# set random seed number for reproducability 
seed = 88;

# setting up simulation parameters
N = 20; # number of cells
m1 = 1; # number of springs per cell
m2 = 4;
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
Discrete_Solution_m1, Discrete_Solution_m2, Continuum_Solution = ComparisonSim(N,m1,m2,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type, 
                                                                            prolif, death, embed, α, βv, γv, event_δt, seed, Av);

indicies = [1,1,5,5,11,11]
num_cols = 2
f1 = TissueGrowth.DiscVSContDensity_plot_all(Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies, num_cols)
save("m_springs_compare.png",f1)

Density_cmap =  :cool #:rainbow1
Density_Range = (2,7)

f2 = TissueGrowth.plotResults2D(Discrete_Solution_m1.u, Discrete_Solution_m1.Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; [1/\text{μm}]", (1.2,1.2), N, m1)
save("Square_infill_m1_springs.png",f2)
f3 = TissueGrowth.plotResults2D(Discrete_Solution_m2.u, Discrete_Solution_m2.Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; [1/\text{μm}]", (1.2,1.2), N, m2)
save("Square_infill_m2_springs.png",f3)
cmap = :jet
xbound = 1.1
ybound = 1.1
Cbar_min = 0
Cbar_max = 10
f2 = DiscVSContShape_plot(Discrete_Solution_m2, m2, Continuum_Solution, xbound, ybound, cmap, Cbar_min, Cbar_max)
