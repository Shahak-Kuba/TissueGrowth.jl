
using TissueGrowth

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

# scaling factor 
Λ = 100000

# setting up simulation parameters
m = 2 # number of springs per cell
R₀ = 282.095  # shape radius μm
D = [0.0075].*Λ
l₀ = 1.0
#kf = 70#93.13 
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 24 # days
δt = 0.01
btypes = ["circle"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "2sigmoid" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0;        β = 0.0;      Ot = 0.0;
event_δt = δt

# 2D simulations 
sol, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 11);

Density_cmap =  :cool #:rainbow1
Stress_cmap = :winter 

geo = 1
diffusivity = 1

Density_Range = (0.05,0.2)
Stress_Range = (-30, 5)

f = TissueGrowth.plotResults2D(sol[geo][diffusivity].u, sol[geo][diffusivity].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[1/μm]}", (280,280), N, m)
f2 = TissueGrowth.plotResults2D_Quadrant(sol[geo][diffusivity].u, sol[geo][diffusivity].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[1/μm]}", (280,280), N, m)
f3 = TissueGrowth.plotStress2D_Quadrant(sol[geo][diffusivity].u, sol[geo][diffusivity].ψ, Stress_cmap, Stress_Range, L"\text{Stress} \; ψ \; \text{[N/μm²]}", (280,280))
save("Experimental_Square_Pore_Nonlinear.png", f)
save("Experimental_Square_Pore_Hookean_Quadrant_f0_resting.png", f2)
save("Experimental_Square_Pore_Hookean_Quadrant_Stress_$l₀.png", f3)



## Cell length histogram plotting
l₀ = 10.0

sol_hookean, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 20, ρ_lim);

# ~ change force law in "CellMechanics.jl"

sol_nonlinear, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 20, ρ_lim);

INDEX = 20;
hookean_cell_lengths = 1 ./ sol_hookean[geo][diffusivity].Density[INDEX].data
nonlinear_cell_lengths = 1 ./ sol_nonlinear[geo][diffusivity].Density[INDEX].data

f4 = TissueGrowth.plotForceLawCompareHistogram(hookean_cell_lengths, nonlinear_cell_lengths, L"\text{Cell length} \text{[μm]}")



















# Compare with regression model from Buenzli et al. 2020

# Regression Model Buenzli et al. 2020 equation (1)

Tb = 28.46 # ± 2.00
v = 2.02 # ± 0.22
t = LinRange(0,28.46,1000)

Ω_estimate = 1 .- (t./Tb).^v

# Normalising our approximated Ω from discrete simulation

Ω = sols2D[1][1].Ω
Ωnorm_Discrete = Ω./Ω[1]
t_sim = sols2D[1][1].t

# Analytic Solution
Ωnorm_Analytic = TissueGrowth.Ω_analytic(Ω[1],N,kf,t)./Ω[1]


f3 = TissueGrowth.plotCompareRegressionBuenzli(Ω_estimate, t, Ωnorm_Analytic, t, Ωnorm_Discrete, t_sim)



