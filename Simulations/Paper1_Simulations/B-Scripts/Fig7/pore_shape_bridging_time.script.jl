
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
N1 = 120 #Int(P*q₀) # number of cells
kf = KF/N1
l_min = 5
l_max = 20

# set random seed number for reproducability 
seed = 99

####### 500μm side length square pores ############
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
Tmax = 28 # days
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
#sol_hookean, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"hookean",dist_type,
                    #prolif, death, embed, α, β, Ot, event_δt, seed, 31);

sol_nonlinear_500, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N1,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

Density_cmap =  :cool #:rainbow1
Stress_cmap = :winter 
Density_Range = (0.05,0.2)
axisTicks = [-600, -400, -200, 0, 200, 400, 600]

f1_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear_500[1].u, sol_nonlinear_500[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (350,350), axisTicks, N1, m, 10)
f2_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear_500[2].u, sol_nonlinear_500[2].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (350,350), axisTicks, N1, m, 10)

save("fig8_square_500.png", f1_nonlinear)
save("fig8_hex_500.png", f2_nonlinear)


####### 750μm side length square pores ############
# setting up simulation parameters
N2 = 150
R₀ = 423.1421876608172  # shape radius μm
Tmax = 48

# 2D simulations
sol_nonlinear_750, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N2,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

f1_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear_750[1].u, sol_nonlinear_750[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (500,500), axisTicks, N2, m, 10)
f2_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear_750[2].u, sol_nonlinear_750[2].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (500,500), axisTicks, N2, m, 10)

save("fig8_square_750.png", f1_nonlinear)
save("fig8_hex_750.png", f2_nonlinear)                    

####### 1000μm side length square pores ############
# setting up simulation parameters
N3 = 204
R₀ = 564.1895835477563  # shape radius μm
Tmax = 48

# 2D simulations
sol_nonlinear_1000, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N3,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);


f1_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear_1000[1].u, sol_nonlinear_1000[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (650,650), axisTicks, N3, m, 10)
f2_nonlinear = TissueGrowth.plotResults2D(sol_nonlinear_1000[2].u, sol_nonlinear_1000[2].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (650,650), axisTicks, N3, m, 10)

save("fig8_square_1000.png", f1_nonlinear)
save("fig8_hex_1000.png", f2_nonlinear)       

## simulations for hex when they share the same initial Density
#ρ₀ = 0.05

#R₀ = 282.095  # shape radius μm
#N_hex_500 = 96 # Int(floor(6*√((2/(3*√3))*π*R₀^2)/20))


# 2D simulations
#sol_nonlinear_hex_ρ_500, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N_hex_500,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,["hex"],"nonlinear",dist_type,
#                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

#R₀ = 423.1421876608172   # shape radius μm
#N_hex_750 = 138 #Int(floor(6*√((2/(3*√3))*π*R₀^2)/20))

# 2D simulations
#sol_nonlinear_hex_ρ_750, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N_hex_750,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,["hex"],"nonlinear",dist_type,
#                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);


#R₀ = 564.1895835477563  # shape radius μm
#N_hex_1000 = 186 #Int(floor(6*√((2/(3*√3))*π*R₀^2)/20))

# 2D simulations
#sol_nonlinear_hex_ρ_1000, 🥔, 🌻 = TissueGrowth.GrowthSimulation(N_hex_1000,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,["hex"],"nonlinear",dist_type,
#                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);




f = TissueGrowth.plotMultiAreaVsTime(sol_nonlinear_1000[1].t,sol_nonlinear_750[1].t,sol_nonlinear_500[1].t,sol_nonlinear_1000[1].Ω,sol_nonlinear_1000[2].Ω,sol_nonlinear_750[1].Ω,sol_nonlinear_750[2].Ω,sol_nonlinear_500[1].Ω,sol_nonlinear_500[2].Ω,N3,N2,N1,kf)
save("fig8_area_compare.png", f)

#∇_500 = ((sol_nonlinear_500[1].Ω[end] - sol_nonlinear_500[1].Ω[1])/Tmax)/N1
#∇_750 = ((sol_nonlinear_750[1].Ω[end] - sol_nonlinear_750[1].Ω[1])/Tmax)/N2
#∇_1000 = ((sol_nonlinear_1000[1].Ω[end] - sol_nonlinear_1000[1].Ω[1])/Tmax)/N3


# Time to bridge based on side length
#Tbₛ = (sₛ, kf, q₀) -> sₛ./(4*kf*q₀)
#Tbₕ = (sₕ, kf, q₀) -> (√3 .* sₕ)./(4*kf*q₀)

#Ω₀ = [250000, 100000, 50000, 20000, 10000]
#Sₛ = sqrt.(Ω₀)
#Sₕ = sqrt.(2/(3√3).*Ω₀)

#ratio = Sₛ ./ Sₕ

#q₀ = 1/20

#Tb_square = Tbₛ(Sₛ, kf, q₀)
#Tb_hex = Tbₕ(Sₛ, kf, q₀)

