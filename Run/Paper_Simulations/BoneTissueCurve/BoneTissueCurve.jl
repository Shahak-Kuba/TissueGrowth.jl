using TissueGrowth

# set random seed number for reproducability 
seed = 99

# setting up simulation parameters
N = 100 # number of cells
m = 1 # number of springs per cell
R₀ = 1.05  # shape radius (does not matter for 1D)
D = [0.5, 0.05, 0.005, 0.0005]
l₀ = 0.02
kf = 0.000316
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "1D"
Tmax = 21.0 # days
δt = 0.01
btypes = ["InvertedBellCurve"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
β = 0.0001;        γ = 0.005;      Ot = 62.5;
event_δt = 0.05

cmap = :spring
cmap2 = :jet
cmap3 = :RdBu_6
geo = 1
crange = (30,80)
crange2 = (0.01, 0.03)
crange3 = (-5, 5)

# NOTE: Make sure to change force law in Model/CellMechanics.jl between simulations

# HOOKEAN SPRINGS SIMULATION

sol_hookean = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                prolif, death, embed, β, γ, Ot, event_δt, seed, 5);

diffusivity = 1

fig_interface_hookean = TissueGrowth.plotResults1D(sol_hookean[geo][diffusivity][1].u, sol_hookean[geo][diffusivity][1].Density, cmap, crange, "Density q", D[diffusivity], kf, m, N)
save("Bone_Hookean_D_0.5.png", fig_interface_hookean)

fig_interface_hookean_stationary_bounds = TissueGrowth.plotStationaryBoundary(sol_hookean[geo][diffusivity][1].u, sol_hookean[geo][diffusivity][1].Vₙ, cmap2, crange2, "velocity [mm day⁻¹]")
save("Bone_Hookean_Stationary_Bounds_D_0.5.png", fig_interface_hookean_stationary_bounds)

fig_hookean_stress_time = TissueGrowth.plotThetaVsTime1D(sol_hookean[geo][diffusivity][1].u, sol_hookean[geo][diffusivity][1].t, 
            sol_hookean[geo][diffusivity][1].ψ, cmap3, crange3, "Stress ψ")
save("Bone_Hookean_Stress_1D.png", fig_hookean_stress_time)

# NONLINEAR SPRINGS SIMULATION

sol_nonlinear = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                prolif, death, embed, β, γ, Ot, event_δt, seed, 5);

diffusivity = 3

fig_interface_nonlinear = TissueGrowth.plotResults1D(sol_nonlinear[geo][diffusivity][1].u, sol_nonlinear[geo][diffusivity][1].Vₙ, cmap, crange, "Density q", D[diffusivity], kf, m, N)
save("Bone_Nonlinear_1D.png", fig_interface_nonlinear)

fig_interface_nonlinear_stationary_bounds = TissueGrowth.plotStationaryBoundary(sol_nonlinear[geo][diffusivity][1].u, sol_nonlinear[geo][diffusivity][1].Vₙ, cmap2, crange2, "velocity [mm day⁻¹]")
save("Bone_Nonlinear_Stationary_Bounds_D_0.005.png", fig_interface_nonlinear_stationary_bounds)

fig_nonlinear_stress_time = TissueGrowth.plotThetaVsTime1D(sol_nonlinear[geo][diffusivity][1].u, sol_nonlinear[geo][diffusivity][1].t, 
            sol_nonlinear[geo][diffusivity][1].ψ, cmap3, crange3, "Stress ψ")
save("Bone_Nonlinear_Stress_D_0.005.png", fig_nonlinear_stress_time)


## HISTOGRAM PLOTTING
timestep = 2
g = TissueGrowth.plotForceLawCompare1D(1 ./ sol_nonlinear[geo][diffusivity][1].Density[timestep][2:end-1], 1 ./ sol_hookean[geo][diffusivity][1].Density[timestep][2:end-1], "Length")
fig_name = "cell_length_$timestep.png"
save(fig_name, g)