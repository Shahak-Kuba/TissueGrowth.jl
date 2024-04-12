using TissueGrowth

# set random seed number for reproducability 
seed = 99

# setting up simulation parameters
N = 120 # number of cells
m = 2 # number of springs per cell
R₀ = 1.05  # shape radius (does not matter for 1D)
D = [0.005]
l₀ = 0.02
kf = 0.001
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 21.0 # days
δt = 0.01
btypes = ["square"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
β = 0.0001;        γ = 0.005;      Ot = 62.5;
event_δt = 0.05

cmap = :spring
cmap2 = :jet
cmap3 = :RdBu_6
geo = 1
diffusivity = 1
crange = (10,40)
crange3 = (0, 35)

sol_hookean = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                prolif, death, embed, β, γ, Ot, event_δt, seed, 6);

fig_interface_hookean = TissueGrowth.plotResults2D(sol_hookean[geo][diffusivity][1].u, sol_hookean[geo][diffusivity][1].Density, cmap2, crange,  "Density q", (1.2,1.2), N, m)
save("Pore_Hookean_Square.png", fig_interface_hookean)

fig_hookean_stress_time = TissueGrowth.plotThetaVsTime(sol_hookean[geo][diffusivity][1].u, sol_hookean[geo][diffusivity][1].t, 
                            sol_hookean[geo][diffusivity][1].ψ, cmap3, crange3, "Stress ψ")
save("Pore_Hookean_Stress_D_0.005.png", fig_hookean_stress_time)


# Make sure to change force law in Model/CellMechanics.jl

sol_nonlinear = TissueGrowth.GrowthSimulation(N,m,R₀,D,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,dist_type,
                prolif, death, embed, β, γ, Ot, event_δt, seed, 8);

fig_interface_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear[geo][diffusivity][1].u, sol_nonlinear[geo][diffusivity][1].Density, cmap, crange, "Density q", (1.2,1.2), N, m)
save("Pore_Nonlinear_Square_Cell_Traj_D_0.005.png", fig_interface_nonlinear)

fig_nonlinear_stress_time = TissueGrowth.plotThetaVsTime(sol_nonlinear[geo][diffusivity][1].u, sol_nonlinear[geo][diffusivity][1].t, 
            sol_nonlinear[geo][diffusivity][1].ψ, cmap3, crange3, "Stress ψ")
save("Pore_Nonlinear_Stress_D_0.005.png", fig_nonlinear_stress_time)


#f = TissueGrowth.plotResults1D(sols1D[geo][diffusivity][1].u, sols1D[geo][diffusivity][1].Density, cmap, crange, "Density q", D[diffusivity], kf, m, N)

timestep = 2
g = TissueGrowth.plotForceLawCompare1D(1 ./ sol_nonlinear[geo][diffusivity][1].Density[timestep][2:end-1], 1 ./ sol_hookean[geo][diffusivity][1].Density[timestep][2:end-1], "Length")
fig_name = "cell_length_$timestep.png"
save(fig_name, g)