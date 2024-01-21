using TissueGrowth

# setting up simulation parameters
N = 77 # number of cells
m = 2 # number of springs per cell
R₀ = 1  # shape radius
D = [0.15] #, 0.075, 0.15, 1] # of cell 
l₀ = 1 # of cell 
kf = 0.01*m # of cell = kf¹
η = 1 # of cell = η¹
Tmax = 25 # days
δt = 0.001
growth_dir = "inward"
btype = "SineWave"
dist_type = "Linear"

sols1D = TissueGrowth.sim1D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type);

cmap = :jet

geo = 1
diffusivity = 1

f = TissueGrowth.plotResults1D(sols1D[diffusivity].u, sols1D[diffusivity].Density, 
                                D[diffusivity], kf,cmap)