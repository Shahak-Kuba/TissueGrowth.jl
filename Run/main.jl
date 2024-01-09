using TissueGrowth


"""
# 1D simulation (periodic boundary)
sols1D = sim1D();

f = plotResults1D(sols1D[1].u, sols1D[1].Density)
f = plotResults1D(sols1D[2].u, sols1D[2].Density)
f = plotResults1D(sols1D[3].u, sols1D[3].Density)
f = plotResults1D(sols1D[4].u, sols1D[4].Density)
SaveData(sols1D, "Alias_2017_redo", "1D_simulations")

f = plotResults1D_Velocity(sols1D[1].u, sols1D[1].Vₙ)
f = plotResults1D_Velocity(sols1D[2].u, sols1D[2].Vₙ)
f = plotResults1D_Velocity(sols1D[3].u, sols1D[3].Vₙ)
f = plotResults1D_Velocity(sols1D[4].u, sols1D[4].Vₙ)

stiffness = 1
plotResults1D_spatial_density(sols1D[stiffness].u, sols1D[stiffness].Vₙ)
#f = plotKapVsVel(sols[1])
"""


# 2D simulations 
sols2D = TissueGrowth.sim2D();
SaveData(sols2D, "Diffusivity_Simulations", "2D_Simulations")
f = plotAreaVStime(sols2D[1])

cmap = :jet

f = TissueGrowth.plotResults2D_Velocity(sols2D[1][1].u, sols2D[1][1].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[2][1].u, sols2D[2][1].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[3][1].u, sols2D[3][1].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[4][1].u, sols2D[4][1].Vₙ, cmap)

f = plotResults2D_Velocity(sols2D[1][2].u, sols2D[1][2].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[2][2].u, sols2D[2][2].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[3][2].u, sols2D[3][2].Vₙ, cmap)

f = plotResults2D_Velocity(sols2D[1][3].u, sols2D[1][3].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[2][3].u, sols2D[2][3].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[3][3].u, sols2D[3][3].Vₙ, cmap)

f = plotResults2D_Velocity(sols2D[1][4].u, sols2D[1][4].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[2][4].u, sols2D[2][4].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[3][4].u, sols2D[3][4].Vₙ, cmap)

f = plotResults2D_Velocity(sols2D[1][5].u, sols2D[1][5].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[2][5].u, sols2D[2][5].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[3][5].u, sols2D[3][5].Vₙ, cmap)

f = plotResults2D_Velocity(sols2D[1][6].u, sols2D[1][6].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[2][6].u, sols2D[2][6].Vₙ, cmap)
f = plotResults2D_Velocity(sols2D[3][6].u, sols2D[3][6].Vₙ, cmap)
# very low diffusivity
f = plotResults2D(sols2D[1][1].u, sols2D[1][1].Density, cmap)
f = plotResults2D(sols2D[1][2].u, sols2D[1][2].Density, cmap)
f = plotResults2D(sols2D[1][3].u, sols2D[1][3].Density, cmap)
f = plotResults2D(sols2D[1][4].u, sols2D[1][4].Density, cmap)
f = plotResults2D(sols2D[1][5].u, sols2D[1][5].Density, cmap)
f = plotResults2D(sols2D[1][6].u, sols2D[1][6].Density, cmap)

# low diffusivity
f = plotResults2D(sols2D[2][1].u, sols2D[2][1].Density, cmap)
f = plotResults2D(sols2D[2][2].u, sols2D[2][2].Density, cmap)
f = plotResults2D(sols2D[2][3].u, sols2D[2][3].Density, cmap)
f = plotResults2D(sols2D[2][4].u, sols2D[2][4].Density, cmap)
f = plotResults2D(sols2D[2][5].u, sols2D[2][5].Density, cmap)
f = plotResults2D(sols2D[2][6].u, sols2D[2][6].Density, cmap)

# mid diffusivity
f = plotResults2D(sols2D[3][1].u, sols2D[3][1].Density, cmap)
f = plotResults2D(sols2D[3][2].u, sols2D[3][2].Density, cmap)
f = plotResults2D(sols2D[3][3].u, sols2D[3][3].Density, cmap)
f = plotResults2D(sols2D[3][4].u, sols2D[3][4].Density, cmap)
f = plotResults2D(sols2D[3][5].u, sols2D[3][5].Density, cmap)
f = plotResults2D(sols2D[3][6].u, sols2D[3][6].Density, cmap)

# high diffusivity
f = plotResults2D(sols2D[4][1].u, sols2D[3][1].Density, cmap)
f = plotResults2D(sols2D[4][2].u, sols2D[3][2].Density, cmap)
f = plotResults2D(sols2D[4][3].u, sols2D[3][3].Density, cmap)
f = plotResults2D(sols2D[4][4].u, sols2D[3][4].Density, cmap)
f = plotResults2D(sols2D[4][5].u, sols2D[3][5].Density, cmap)
f = plotResults2D(sols2D[4][6].u, sols2D[3][6].Density, cmap)

# area loss simulations (take around 1 hr)
sols2d_δt = sim2D_δt();
plotAreaVStime_δt(sols2d_δt[1][1], sols2d_δt[1][2], sols2d_δt[1][3], sols2d_δt[1][4], "triangle")
plotAreaVStime_δt(sols2d_δt[2][1], sols2d_δt[2][2], sols2d_δt[2][3], sols2d_δt[2][4], "square")
plotAreaVStime_δt(sols2d_δt[3][1], sols2d_δt[3][2], sols2d_δt[3][3], sols2d_δt[3][4], "hex")
plotAreaVStime_δt(sols2d_δt[4][1], sols2d_δt[4][2], sols2d_δt[4][3], sols2d_δt[4][4], "star")
plotAreaVStime_δt(sols2d_δt[5][1], sols2d_δt[5][2], sols2d_δt[5][3], sols2d_δt[5][4], "cross")
SaveData(sols2d_δt, "AreaLoss_Simulations", "2D_Simulations")

cmap = :turbo
geo = 1
dt_sim = 4
f = plotResults2D(sols2d_δt[geo][dt_sim].u, sols2d_δt[geo][dt_sim].Density, cmap)