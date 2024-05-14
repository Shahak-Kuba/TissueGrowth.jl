
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
m = 1 # number of springs per cell
m2 = 2
R₀ = 282.095  # shape radius μm
D = 0.00
kₛ = 7.5
Kₛ = 200
l₀ = 20.0
L₀ = ((l_max - l_min)/((ks/Ks)*((l_max^2 - l_min^2)/2 + a_hookean*(l_min - l_max)) - log(l_min/l_max)))
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 24 # days
δt = 0.01
btypes = ["square"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0;        β = 0.0;      Ot = 0.0;
event_δt = δt

# 2D simulations
sol_hookean, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"hookean",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 5);

sol_nonlinear, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 5);

# 2D simulations (stress)
sol_hookean2, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m2,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"hookean",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

sol_nonlinear2, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m2,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

Density_cmap =  :cool #:rainbow1
Stress_cmap = :winter 

geo = 1
Stress_Range_Hookean = (-2, 2)
Stress_Range_Nonlinear = (48, 52)

f_length_Hookean = TissueGrowth.plotForceLawCompareStairs(sol_hookean[geo].Density)
f2_hookean = TissueGrowth.plotStress2D_Quadrant(sol_hookean2[geo].u, sol_hookean2[geo].ψ, Stress_cmap, Stress_Range_Hookean, L"σ/E \; \text{[-]}", (300,300))

f_length_Nonlinear = TissueGrowth.plotForceLawCompareStairs(sol_nonlinear[geo].Density)
f2_nonlinear= TissueGrowth.plotStress2D_Quadrant(sol_nonlinear2[geo].u, sol_nonlinear2[geo].ψ, Stress_cmap, Stress_Range_Nonlinear, L"σ/E \; \text{[-]}", (300,300))

save("Fig7_High_Hookean_Cell_Length.png", f_length_Hookean)
save("Fig7_High_Hookean_Stress.png", f2_hookean)

save("Fig7_High_Nonlinear_Cell_Length.png", f_length_Nonlinear)
save("Fig7_High_Nonlinear_Stress.png", f2_nonlinear)