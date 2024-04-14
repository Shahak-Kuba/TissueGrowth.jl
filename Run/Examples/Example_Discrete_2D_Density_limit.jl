using TissueGrowth

# set random seed number for reproducability 
seed = 99

# setting up simulation parameters
N = 120 # number of cells
m = 1 # number of springs per cell
R₀ = 1.05  # shape radius
D = [0.01]
l₀ = 1.0
kf = 0.0013
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 22.0 # days
δt = 0.01
btypes = ["square"]#, "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        β = 0.001;      Ot = 62.5;
event_δt = δt

# Cell density limit
q_lim = 40

# 2D simulations 
sol, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 11, q_lim);

Density_cmap = :jet
crange = (10, 60)
ψ_cmap = :balance

geo = 1
diffusivity = 1

Density_Range = (10,30)
ψ_Range = (-20,0)
f = TissueGrowth.plotResults2D(sol[geo][diffusivity].u, sol[geo][diffusivity].Density, Density_cmap, crange,  "Density q", (1.2,1.2), N, m)
f = TissueGrowth.plotResults2D(sol[diffusivity][geo].u, sol[diffusivity][geo].t, sol[diffusivity][geo].ψ, ψ_cmap, ψ_Range, "Stress ψ", D[diffusivity], kf)
