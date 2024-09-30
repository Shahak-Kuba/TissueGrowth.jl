using TissueGrowth
using Makie

# set random seed number for reproducability 
seed = 42

# See parameter approximation document
# Calculating kf
KF = 8784.2;
Tb = 28.46
l = 500;
Ω₀ = l^2
P = l*4
q₀ = 1/20; 
N = 100 #Int(P*q₀) # number of cells
kf = KF/N
l_min = 10
l_max = 20


# setting up simulation parameters
m = 6 # number of springs per cell
R₀ = 282.095  # shape radius μm
D = 0.00
Kₛ = 12
L₀ = 10.0
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 30 # days
δt = 0.01
btypes = ["square"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m
restoring_force = "nonlinear"

## Cell Behaviours
prolif = true; death = false; embed = true;
α = 0.2;        β = 0.0;      Ot = 0.0001;
event_δt = δt

# 2D simulations 
sols2D, embedded_cells, embed_cell_count = TissueGrowth.GrowthSimulation(N,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,restoring_force,dist_type,
                                            prolif, death, embed, α, β, Ot, event_δt, seed, 200);

geo = 1
diffusivity = 1

Density_cmap =  :cool #:rainbow1
Density_Range = (0.01,0.4)

f2 = TissueGrowth.plotResults2D(sols2D[1].u, sols2D[1].Density, Density_cmap, Density_Range,  "Density q μm⁻¹", (300,300), N, m, 20)

f3 = TissueGrowth.plotResults2D_embedded(sols2D[1].u, sols2D[1].Density, Density_cmap, Density_Range, "q [1/μm]", D, kf, (300,300), embedded_cells, true)
save("Cell_culturing_Interface_Plot.png",f3)

