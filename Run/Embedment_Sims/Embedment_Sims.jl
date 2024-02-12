using TissueGrowth

# set random seed number for reproducability 
seed = 211

# setting up simulation parameters
N = 60 # number of cells
m = 3 # number of springs per cell
R₀ = 1.05  # shape radius
D = [0.005]
l₀ = 1.0
kf = 0.002
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
Tmax = 39.0 # days
δt = 0.01
btypes = ["square"]#, "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = true;
α = 0.0001;        β = 0.001;      Ot = 30;
event_δt = δt

# 2D simulations 
sols2D, embedded_cells, embed_cell_count = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 100);

geo = 1
diffusivity = 1

f = TissueGrowth.plotOtValueVsTime(sols2D[diffusivity][geo].t, sols2D[diffusivity][geo].Ω, embed_cell_count[1], Ot/m)

# plotting interface
Density_Range = (20,40)
Density_cmap = :jet
multiple_Interfaces = false
f_interface = TissueGrowth.plotResults2D(sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].Density, Density_cmap, Density_Range, 
                                        "Density ρ", D[diffusivity], kf, embedded_cells, multiple_Interfaces)