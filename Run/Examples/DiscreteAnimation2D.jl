using TissueGrowth

# set random seed number for reproducability 
seed = 88

# setting up simulation parameters
N = 60 # number of cells
m = 2 # number of springs per cell
R₀ = 1.0  # shape radius
D = [0.01]
l₀ = 1.0
kf = 0.002
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
Tmax = 26.0 # days
δt = 0.01
btypes = ["square"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        β = 0.001;      γ = 0.01;
event_δt = δt

NumSaveTimePoints = 100;

# 2D simulations 
sols2D = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,
            prolif, death, embed, α, β, γ, event_δt, seed, NumSaveTimePoints);

Density_cmap = :jet

geo = 1
diffusivity = 1

Density_Range = (10,30)
filename = "squarePore.gif"

# making Animation
TissueGrowth.animateResults2D(sols2D[diffusivity][geo].t, sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].Density, Density_cmap, Density_Range, "Density ρ", D[diffusivity], kf, filename)