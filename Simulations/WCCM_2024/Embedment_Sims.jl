using TissueGrowth

# set random seed number for reproducability 
seed = 77

# See parameter approximation document
# Calculating kf
KF = 8784.2;
Tb = 28.46
l = 500;
Ω₀ = l^2
P = l*4
q₀ = 1/20; 
N = Int(P*q₀) # number of cells
kf = KF/N
l_min = 5
l_max = 20


# setting up simulation parameters
m = 3 # number of springs per cell
R₀ = 100 #282.095  # shape radius μm
D = 0.00
kₛ = 7.5
Kₛ = 60
l₀ = 10.0
L₀ = ((l_max - l_min)/((kₛ/Kₛ)*((l_max^2 - l_min^2)/2 + l₀*(l_min - l_max)) - log(l_min/l_max)))
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 15 # days
δt = 0.01
btypes = ["PerturbedCircle"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m

## Cell Behaviours
prolif = false; death = false; embed = true;
α = 0.0;        β = 0.0;      Ot = 0.000625;
event_δt = δt

# 2D simulations 
sols2D, embedded_cells, embed_cell_count = TissueGrowth.GrowthSimulation(N,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                                            prolif, death, embed, α, β, Ot, event_δt, seed, 200);

geo = 1
diffusivity = 1

Density_cmap =  :cool #:rainbow1
Density_Range = (0.05,0.15)
f2 = TissueGrowth.plotResults2D(sols2D[1].u, sols2D[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; [\text{μm^{-1}}]", (250,250), N, m, 20)

f3 = TissueGrowth.plotResults2D_embedded(sols2D[1].u, sols2D[1].Density, Density_cmap, Density_Range, L"\text{Density} \; q \; [\text{μm^{-1}}]", D, kf, (250,250), embedded_cells, true)

f = TissueGrowth.plotOtValueVsTime(sols2D[1].t, sols2D[1].Ω, embed_cell_count[1], Ot/m, m)

# plotting interface
Density_Range = (20,40)
Density_cmap = :jet
multiple_Interfaces = false
f_interface = TissueGrowth.plotResults2D(sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].Density, Density_cmap, Density_Range, 
                                        "Density ρ", D[diffusivity], kf, embedded_cells, multiple_Interfaces)