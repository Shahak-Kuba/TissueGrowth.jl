# Simulations based off Lanaro et al. 2021 Experimental Investigation
# Square scaffold side length = 500 μm → Perimeter = 2000 μm, Area = 250000 μm²
# Pore radius (circular pore) = √(Area / π) ≈ 282.09 μm

# MC3T3-E1 cell diameter ≈ 20 to 50 μm,  10.1089/ten.tea.2011.0545, Gibon et al. 2012
# Assuming 20 μm cell diameter ⟹ 2000 μm / 20 μm = 100 cells along the scaffold at the initial time.

# 500 μm square scaffold can fill within 28 days (see Fig 2 Lanaro et al. 2021)
# Assuming kf is constant for all cells, kf = (Pore Area / Time to close / Number of cells) = (2500 μm² / 28 days / 100 cells) ≈ 0.89 μm²/day/cell

# FOR SIMULATION
# To maintain cells under tension, assume closed pore radius = 50 μm (to avoid high densities and numerical instabilities) → void area ≈ 7853.98
# Given constant number of cells, the ending cell length = 2π(close radius)/Number of cells = π
# Assume closing time < 28 days ≈ 26 days
# kf = (Filled Area) / Time to close / Number of cells = ((250000 - 7853.98)) / 26 / 100 ≈ 93.13
# initial density = 1 / 20 μm = 0.05 , final density = 1 / π ≈ 0.32
using TissueGrowth

# set random seed number for reproducability 
seed = 99

# scaling factor to simulate μm instead of mm
Λ = 1000

# setting up simulation parameters
N = 100 # number of cells
m = 4 # number of springs per cell
R₀ = 282.09  # shape radius μm
D = [0.5].*Λ
l₀ = 3.14
kf = 93.13 
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
Tmax = 26.0 # days
δt = 0.01
btypes = ["square"] #, "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        β = 0.001;      Ot = 62.5;
event_δt = δt

# 2D simulations 
sols2D, 🥔, 🌻 = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 11);

Density_cmap = :jet

geo = 1
diffusivity = 1

Density_Range = (0.05,0.32)

f = TissueGrowth.plotResults2D(sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].Density, Density_cmap, Density_Range, "Density ρ", D[diffusivity], kf, (500,500))