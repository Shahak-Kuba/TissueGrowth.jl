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


# setting up simulation parameters
m = 1 # number of springs per cell
R₀ = 282.095  # shape radius μm
D = 0.00
kₛ = 500
l₀ = 10
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 26 # days
δt = 0.01
btypes = ["square"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m

## Cell Behaviours
prolif = false; death = false; embed = false;
β = 0.0;        γ = 0.0;      Ot = 0.0;
event_δt = δt

cmap = :cool
cmap2 = :jet
cmap3 = :RdBu_6
geo = 1
diffusivity = 1
crange = (0.05,0.15)
crange3 = (0, 35)

sol_hookean = TissueGrowth.GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                prolif, death, embed, β, γ, Ot, event_δt, seed, 8);

fig_interface_hookean = TissueGrowth.plotResults2D(sol_hookean[1][1].u, sol_hookean[1][1].Density, cmap, crange,  "Density q", (280,280), N, m, 8)
save("Pore_Hookean_Square.png", fig_interface_hookean)

fig_hookean_stress_time = TissueGrowth.plotThetaVsTime(sol_hookean[geo][diffusivity][1].u, sol_hookean[geo][diffusivity][1].t, 
                            sol_hookean[geo][diffusivity][1].ψ, cmap3, crange3, "Stress ψ")
save("Pore_Hookean_Stress_D_0.005.png", fig_hookean_stress_time)


# Make sure to change force law in Model/CellMechanics.jl

sol_nonlinear = TissueGrowth.GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                prolif, death, embed, β, γ, Ot, event_δt, seed, 5);

fig_interface_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear[1][1].u, sol_nonlinear[1][1].Density, cmap, crange, "Density q", (280,280), N, m, 5)
save("Pore_Nonlinear_Square_Cell_Traj_D_0.005.png", fig_interface_nonlinear)

fig_interface_nonlinear_Quadrant = TissueGrowth.plotResults2D_Quadrant(sol_nonlinear[geo][diffusivity][1].u, sol_nonlinear[geo][diffusivity][1].Density, cmap, crange, "Density q", (1.2,1.2), N, m)

fig_nonlinear_stress_time = TissueGrowth.plotThetaVsTime(sol_nonlinear[geo][diffusivity][1].u, sol_nonlinear[geo][diffusivity][1].t, 
            sol_nonlinear[geo][diffusivity][1].ψ, cmap3, crange3, "Stress ψ")
save("Pore_Nonlinear_Stress_D_0.005.png", fig_nonlinear_stress_time)


idx = 3
t = round(sol_nonlinear[1][1].t[idx], digits=1)
h = plotForceLawCompareStairs(sol_nonlinear[1][1].Density)
save("Cell_length_hist_t_$t.png", h)
