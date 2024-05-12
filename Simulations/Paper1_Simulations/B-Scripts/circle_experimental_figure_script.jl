
using TissueGrowth
using Makie

# See parameter approximation document
# Calculating kf
KF = 8784.2;
Tb = 28.46
l = 500;
Ω₀ = l^2
P = l*4
V₀ = abs.(TissueGrowth.V(KF, Ω₀, 0))
q₀ = 1/20;
N = Int(P*q₀) # number of cells
kf = KF/N

# set random seed number for reproducability 
seed = 99


# setting up simulation parameters
m = 2 # number of springs per cell
R₀ = 282.095  # shape radius μm
D = 0.00
kₛ = 7.5
l₀ = 10 
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 24 # days
δt = 0.01
btypes = ["circle"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0;        β = 0.0;      Ot = 0.0;
event_δt = δt

# 2D simulations
sol_hookean, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"hookean",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

sol_nonlinear, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

Density_cmap =  :cool #:rainbow1
Stress_cmap = :winter 

geo = 1

Density_Range = (0.05,0.2)
Stress_Range = (-50, 10)

f1_hookean = TissueGrowth.plotResults2D_Quadrant(sol_hookean[geo].u, sol_hookean[geo].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[1/μm]}", (280,280), N, m, 10)
f2_hookean = TissueGrowth.plotStress2D_Quadrant(sol_hookean[geo].u, sol_hookean[geo].ψ, Stress_cmap, Stress_Range, L"\text{Stress} \; σ \; \text{[N/μm²]}", (280,280))

f1_nonlinear = TissueGrowth.plotResults2D_Quadrant(sol_nonlinear[geo].u, sol_nonlinear[geo].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[1/μm]}", (280,280), N, m, 10)
f2_nonlinear = TissueGrowth.plotStress2D_Quadrant(sol_nonlinear[geo].u, sol_nonlinear[geo].ψ, Stress_cmap, Stress_Range, L"\text{Stress} \; σ \; \text{[N/μm²]}", (280,280))

save("circle_hookean_cell_traj.png", f1_hookean)
save("circle_hookean_stress.png", f2_hookean)

save("circle_nonlinear_cell_traj.png", f1_nonlinear)
save("circle_nonlinear_stress.png", f2_nonlinear)



