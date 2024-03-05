using TissueGrowth

# set random seed number for reproducability 
seed = 99

# setting up simulation parameters
N = 120 # number of cells
m = 1 # number of springs per cell
R₀ = 1.05  # shape radius
D = [0.001, 0.05, 1]
l₀ = 1.0
kf = 0.001
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
Tmax = 24.0 # days
δt = 0.01
btypes = ["square", "hex"]#, "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        β = 0.001;      Ot = 62.5;
event_δt = δt

# 2D simulations 
sols2D, 🥔, 🌻 = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 11);


# Plotting
axislims = (1.2,1.2)
cmap = :jet
CRange = (15, 40)
TissueGrowth.plotMultiSimResults2D(sols2D, axislims, cmap, CRange)
