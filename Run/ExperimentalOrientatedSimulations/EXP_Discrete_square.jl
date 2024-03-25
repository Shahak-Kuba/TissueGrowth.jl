# See parameter approximation document


using TissueGrowth

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

# scaling factor to simulate μm instead of mm
Λ = 10000

# setting up simulation parameters
m = 1 # number of springs per cell
R₀ = 282.095  # shape radius μm
D = [0.05].*Λ
l₀ = 3.14
#kf = 70#93.13 
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
Tmax = 28.4 # days
δt = 0.01
btypes = ["circle"] #, "triangle", "square", "hex", "star","cross"] #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0;        β = 0.0;      Ot = 0.0;
event_δt = δt

# 2D simulations 
sols2D, 🥔, 🌻 = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,
            prolif, death, embed, α, β, Ot, event_δt, seed, 140);

Density_cmap = :jet
Stress_cmap = :viridis

geo = 1
diffusivity = 1

Density_Range = (0.05,0.32)
Stress_Range = (-20, 20)

f = TissueGrowth.plotResults2D(sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].Density, Density_cmap, Density_Range, "Density ρ", D[diffusivity], kf, (280,280))
f2 = TissueGrowth.plotThetaVsTime(sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].t, sols2D[diffusivity][geo].ψ, Stress_cmap, Stress_Range, "Stress ψ", D, kf)

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



