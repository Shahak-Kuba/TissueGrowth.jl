using TissueGrowth
include("ComparisonPlottingFncs.jl")

# Shared variables
R₀ = 1.05  
D = 0.005
kf = 0.001
growth_dir = "inward"
Tmax = 26.0 # days
btype = "square" #Options: ["circle", "triangle", "square", "hex", "star","cross"]

# Discrete Simulation
# set random seed number for reproducability 
seed = 88

# setting up simulation parameters
N = 96 # number of cells
m1 = 1 # number of springs per cell
l₀ = 1.0
η = 1.0 
δt = 0.01
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        βv = 0.001;      γv = 0.01;
event_δt = δt

# 2D simulations 
sols2D_m1 = TissueGrowth.sim2D(N,m1,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,[btype],dist_type,
            prolif, death, embed, α, βv, γv, event_δt, seed, 11);

m2 = 4 # number of springs per cell
sols2D_m2 = TissueGrowth.sim2D(N,m2,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,[btype],dist_type,
            prolif, death, embed, α, βv, γv, event_δt, seed, 11);


# Continuum Simulation
Av = 0.00;
ρ₀ = sols2D_m1[1][1].Density[1][1];

θ_cont,R_cont,ρ_cont = TissueGrowth.SolveContinuumLim_Polar(D,kf,Av,ρ₀,Tmax,R₀,btype,growth_dir);



# plotting
Discrete_Solution_m1 = sols2D_m1[1][1];
Discrete_Solution_m2 = sols2D_m2[1][1];
Continuum_Solution = (θ_cont,R_cont,ρ_cont);
indicies = [1,3,5,7,9,11]
num_cols = 2
DiscVSContDensity_plot_all(Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies, num_cols)