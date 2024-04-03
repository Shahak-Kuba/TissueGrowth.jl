using TissueGrowth

# set random seed number for reproducability 
seed = 99

# setting up simulation parameters
N = 75 # number of cells
m = 1 # number of springs per cell
R₀ = 1.05  # shape radius
D = [0.005]
l₀ = 1.0
kf = 0.00031
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "1D"
Tmax = 30.0 # days
δt = 0.01
btypes = ["InvertedBellCurve"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
β = 0.0001;        γ = 0.005;      Ot = 62.5;
event_δt = 0.05


sols1D = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                prolif, death, embed, β, γ, Ot, event_δt, seed, 30);

cmap = :spring
geo = 1
diffusivity = 1
crange = (30,60)

#f = TissueGrowth.plotResults1D(sols1D[geo][diffusivity][1].u, sols1D[geo][diffusivity][1].Density, 
#                                D[diffusivity], kf,cmap, 12, 7)

f = TissueGrowth.plotResults1D(sols1D[geo][diffusivity][1].u, sols1D[geo][diffusivity][1].Density, cmap, crange, "Density q", D[diffusivity], kf)
f = TissueGrowth.plotResults1D(sols1D[geo][diffusivity][1].u, sols1D[geo][diffusivity][1].Density, cmap, crange, "Density q", D[diffusivity], kf, m, N)
