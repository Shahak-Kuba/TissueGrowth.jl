using TissueGrowth
include("ComparisonPlottingFncs.jl")

# Shared variables
R₀ = 1.05  
D = 0.005
kf = 0.0005
growth_dir = "inward"
Tmax = 26.0 # days
btype = "square" #Options: ["circle", "triangle", "square", "hex", "star","cross"]

# Discrete Simulation
# set random seed number for reproducability 
seed = 88

# setting up simulation parameters
N = 180 # number of cells
m = 1 # number of springs per cell
l₀ = 1.0
η = 1.0 
δt = 0.01
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        β = 0.001;      γ = 0.01;
event_δt = δt

# 2D simulations 
sols2D = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,[btype],dist_type,
            prolif, death, embed, α, β, γ, event_δt, seed, 11);


# Continuum Simulation
A = 0.00;
ρ₀ = sols2D[1][1].Density[1][1];

θ_cont,R_cont,ρ_cont = SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀,btype,growth_dir);



# plotting
Discrete_Solution = sols2D[1][1];
Continuum_Solution = (θ_cont,R_cont,ρ_cont);
index = 9
DiscVSContDensity_plot(Discrete_Solution, Continuum_Solution, index)