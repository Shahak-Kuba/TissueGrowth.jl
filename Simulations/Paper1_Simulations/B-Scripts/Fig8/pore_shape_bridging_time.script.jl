
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
N = 120 #Int(P*q₀) # number of cells
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
Kₛ = 150
l₀ = 10.0
L₀ = ((l_max - l_min)/((kₛ/Kₛ)*((l_max^2 - l_min^2)/2 + l₀*(l_min - l_max)) - log(l_min/l_max)))
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 24 # days
δt = 0.01
btypes = ["square", "hex"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
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


Density_Range = (0.05,0.2)

f1_hookean = TissueGrowth.plotResults2D(sol_hookean[1].u, sol_hookean[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (280,280), N, m, 10)
f2_hookean = TissueGrowth.plotResults2D(sol_hookean[2].u, sol_hookean[2].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (320,320), N, m, 10)

f1_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear[1].u, sol_nonlinear[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (280,280), N, m, 10)
f2_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear[2].u, sol_nonlinear[2].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (320,320), N, m, 10)

save("fig8_square_hookean_cell.png", f1_hookean)
save("fig8_hex_hookean_cell.png", f2_hookean)

save("fig8_square_nonlinear_cell.png", f1_nonlinear)
save("fig8_hex_nonlinear_cell.png", f2_nonlinear)


f = TissueGrowth.plotMultiAreaVsTime(sol_hookean[1].Ω,sol_hookean[1].t,sol_hookean[2].Ω,sol_hookean[1].t,sol_nonlinear[1].Ω,sol_nonlinear[1].t,sol_nonlinear[2].Ω,sol_nonlinear[1].t,N,kf)
save("fig8_area_compare.png", f)


# Time to bridge based on side length
Tbₛ = (sₛ, kf, q₀) -> sₛ./(4*kf*q₀)
Tbₕ = (sₕ, kf, q₀) -> (√3 .* sₕ)./(4*kf*q₀)

Ω₀ = [250000, 100000, 50000, 20000, 10000]
Sₛ = sqrt.(Ω₀)
Sₕ = sqrt.(2/(3√3).*Ω₀)

ratio = Sₛ ./ Sₕ

q₀ = 1/20

Tb_square = Tbₛ(Sₛ, kf, q₀)
Tb_hex = Tbₕ(Sₛ, kf, q₀)

