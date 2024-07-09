using TissueGrowth

# set random seed number for reproducability 
seed = 10

# See parameter approximation document
# Calculating kf
KF = 8784.2;
Tb = 28.46
l = 500;
Ω₀ = l^2
P = l*4
q₀ = 1/20; 
N = 50#Int(P*q₀) # number of cells
kf = KF/N
l_min = 5
l_max = 20


# setting up simulation parameters
m = 6 # number of springs per cell
R₀ = 156.4264450801666  # shape radius μm (To match other pore)
D = 0.00
kₛ = 1
Kₛ = 15
l₀ = 10.0
L₀ = ((l_max - l_min)/((kₛ/Kₛ)*((l_max^2 - l_min^2)/2 + l₀*(l_min - l_max)) - log(l_min/l_max)))
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 8 # days
δt = 0.01
btypes = ["square"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0;        β = 0.0;      Ot = 0.003;
event_δt = δt

# 2D simulations 
sols2D, embedded_cells, embed_cell_count = TissueGrowth.GrowthSimulation(N,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                                            prolif, death, embed, α, β, Ot, event_δt, seed, 200);

geo = 1
diffusivity = 1

Density_cmap =  :cool #:rainbow1
Density_Range = (0.02,0.1)
filename = "WCCM_2024_squarePore.gif"
#f = TissueGrowth.plotOtValueVsTime(sols2D[1].t, sols2D[1].Ω, embed_cell_count[1], Ot/m, m)
#save("WCCM_2024_Square_Ot_Plot.png",f)

f2 = TissueGrowth.plotResults2D(sols2D[1].u, sols2D[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; [\text{μm^{-1}}]", (200,200), N, m, 20)

#f3 = TissueGrowth.plotResults2D_embedded(sols2D[1].u, sols2D[1].Density, Density_cmap, Density_Range, L"\text{Density} \; q \; [\text{μm^{-1}}]", D, kf, (200,200), embedded_cells, true)
#save("WCCM_2024_Square_Interface_Plot.png",f3)

TissueGrowth.animateResults2D(sols2D[1].t, sols2D[1].u, sols2D[1].Density, Density_cmap, Density_Range, "Cell density q", D[diffusivity], kf, filename)
