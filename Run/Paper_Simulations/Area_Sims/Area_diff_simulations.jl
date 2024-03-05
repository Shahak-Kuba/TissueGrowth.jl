using TissueGrowth
using DifferentialEquations

EXP_VALUES = true

# set random seed number for reproducability 
seed = 99

if !EXP_VALUES
    # setting up simulation parameters
    N = 100 # number of cells
    m = 1 # number of springs per cell
    R₀ = 1.05  # shape radius
    D = [0.05]
    l₀ = 1.0
    kf = 0.001
    η = 1.0 
    growth_dir = "inward" # Options: "inward", "outward"
    Tmax = 26.0 # days
    δt_array = [0.01, 0.005, 0.001, 0.0001]
    btypes = ["square"]#, "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
    dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
else
    # scaling factor to simulate μm instead of mm
    Λ = 10000

    # setting up simulation parameters
    N = 100 # number of cells
    m = 1 # number of springs per cell
    R₀ = 282.09  # shape radius μm
    D = [0.05].*Λ
    l₀ = 3.14
    kf = 93.13 
    η = 1.0 
    growth_dir = "inward" # Options: "inward", "outward"
    Tmax = 26.0 # days
    δt = 0.01
    btypes = ["square"] #, "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
    dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

end

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        β = 0.001;      Ot = 62.5;
event_δt = δt_array[1]

# 2D simulations 
sols2D_δt_1, 🥔, 🌻 = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt_array[1],btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 500, DifferentialEquations.Euler());

Ω₁ = sols2D_δt_1[1][1].Ω;
t₁ = sols2D_δt_1[1][1].t;

sols2D_δt_2, 🥔, 🌻 = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt_array[2],btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 500, DifferentialEquations.Euler());

Ω₂ = sols2D_δt_2[1][1].Ω;
t₂ = sols2D_δt_2[1][1].t;

sols2D_δt_3, 🥔, 🌻 = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt_array[3],btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 500, DifferentialEquations.Euler());

Ω₃ = sols2D_δt_3[1][1].Ω;
t₃ = sols2D_δt_3[1][1].t;

sols2D_δt_4, 🥔, 🌻 = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt_array[4],btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 500, DifferentialEquations.Euler());

Ω₄ = sols2D_δt_4[1][1].Ω;
t₄ = sols2D_δt_4[1][1].t;


TissueGrowth.plotδtAreaResults(Ω₁,t₁,Ω₂,t₂,Ω₄,t₄,N,kf)


Density_cmap = :jet

geo = 1
diffusivity = 1

Density_Range = (0.05,0.32)

f = TissueGrowth.plotResults2D(sols2D_δt_1[diffusivity][geo].u, sols2D_δt_1[diffusivity][geo].Density, Density_cmap, Density_Range, "Density ρ", D[diffusivity], kf, (280,280))