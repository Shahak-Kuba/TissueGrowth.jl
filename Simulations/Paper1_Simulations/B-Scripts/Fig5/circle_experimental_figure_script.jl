
using TissueGrowth
using Makie

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

# set random seed number for reproducability 
seed = 99


# setting up simulation parameters
m = 2 # number of springs per cell
R₀ = 282.095  # shape radius μm
D = 0.00
kₛ = 7.5
Kₛ = kₛ / 0.2^2
l₀ = 10.0
L₀ = ((l_max - l_min)/((kₛ/Kₛ)*((l_max^2 - l_min^2)/2 + l₀*(l_min - l_max)) - log(l_min/l_max)))
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

sol_nonlinear, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

Density_cmap =  :cool #:rainbow1
Stress_cmap = :winter 

geo = 1

Density_Range = (0.05,0.1)
Stress_Range_Hookean = (-2, 2)
Stress_Range_Nonlinear = (-2, 2)

f1_hookean = TissueGrowth.plotResults2D_Quadrant(sol_hookean[geo].u, sol_hookean[geo].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[1/μm]}", (300,300), N, m, 10)
#f2_hookean = TissueGrowth.plotStress2D_Quadrant(sol_hookean[geo].u, sol_hookean[geo].ψ, Stress_cmap, Stress_Range_Hookean, L"σ/E \; \text{[-]}", (300,300))

f1_nonlinear = TissueGrowth.plotResults2D_Quadrant(sol_nonlinear[geo].u, sol_nonlinear[geo].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[1/μm]}", (300,300), N, m, 10)
#f2_nonlinear = TissueGrowth.plotStress2D_Quadrant(sol_nonlinear[geo].u, sol_nonlinear[geo].ψ, Stress_cmap, Stress_Range_Nonlinear, L"σ/E \; \text{[-]}", (300,300))

save("circle_hookean_cell_traj.png", f1_hookean)
#save("circle_hookean_stress.png", f2_hookean)

save("circle_nonlinear_cell_traj.png", f1_nonlinear)
#save("circle_nonlinear_stress.png", f2_nonlinear)



