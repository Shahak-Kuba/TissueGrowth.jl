using TissueGrowth

# set random seed number for reproducability 
seed = 99

# setting up simulation parameters
N = 180 # number of cells
m = 1 # number of springs per cell
R₀ = 1.05  # shape radius
D = [0.01]
l₀ = 1.0
kf = 0.001
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
Tmax = 20.0 # days
δt = 0.01
btypes = ["square"]#, "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = true;
α = 0.0001;        β = 0.001;      Ot = 62.5;
event_δt = δt

# 2D simulations 
sols2D, 🥔, 🌻 = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 11);

Density_cmap = :jet
ψ_cmap = :balance

geo = 1
diffusivity = 1

Density_Range = (10,30)
ψ_Range = (-20,0)
f = TissueGrowth.plotResults2D(sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].Density, Density_cmap, Density_Range, "Density ρ", D[diffusivity], kf)
f = TissueGrowth.plotResults2D(sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].t, sols2D[diffusivity][geo].ψ, ψ_cmap, ψ_Range, "Stress ψ", D[diffusivity], kf)
