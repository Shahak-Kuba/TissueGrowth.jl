using TissueGrowth

# set random seed number for reproducability 
seed = 99

# setting up simulation parameters
N = 120 # number of cells
m = 2 # number of springs per cell
R₀ = 1.2694265629824517 #Almie-Hex: 1.3640876152390462 #Almie-Square: 1.2694265629824517  # shape radius (NOTE: Almie was based on perimeter not same initial area)
D = [0.001, 0.005, 1]
l₀ = 1.0
kf = 0.0012
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 26.0 # days
δt = 0.01
btypes = ["square", "hex"]#, "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 100
ρ_lim = q_lim * m

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        β = 0.001;      Ot = 62.5;
event_δt = δt

# 2D simulations 
sols2D, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 8, ρ_lim);


# Plotting
axislims = (1.6,1.6)
cmap = :jet
CRange = (0.01, 0.03)
shape_comp_fig = TissueGrowth.plotMultiSimResults2D(sols2D, axislims, cmap, CRange)
save("Shape_Compare_Square.png", shape_comp_fig)

Ω₁ = sols2D[2][1].Ω
t₁ = sols2D[2][1].t
Ω₂ = sols2D[2][1].Ω
t₂ = sols2D[2][1].t
TissueGrowth.plotMultiAreaVsTime(Ω₁,t₁,Ω₂,t₂,N,kf)

