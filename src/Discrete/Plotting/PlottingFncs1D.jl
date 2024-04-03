function findMinMax(var)
    Min = 0
    Max = 0
    for ii in eachindex(var)
        if ii == 1
            Min = minimum(var[ii])
            Max = maximum(var[ii])
        else
            if minimum(var[ii]) < Min
                Min = minimum(var[ii])
            end
            if maximum(var[ii]) > Max
                Max = maximum(var[ii])
            end
        end
    end
    return (Min, Max)
end

function plotResults1D(u, var, cmap, crange, cbarlabel, D, kf)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(0, 1.5, -0.1, 1.1), aspect=DataAspect(), 
              xlabel="x", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="y", ylabelsize = txtSize, yticklabelsize = tickSize,
              title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    for i in eachindex(u)
        plotInterface1D!(gaxmain, u, var, cmap, CRange, i, 5)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotResults1D(u, var, cmap, crange, cbarlabel, D, kf, m, N)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-0.01, 1.5, -0.1, 1.1), aspect=DataAspect(), 
              xlabel="x", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="y", ylabelsize = txtSize, yticklabelsize = tickSize,
              title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    for i in eachindex(u)
        plotInterface1D!(gaxmain, u, var, cmap, CRange, i, 5)
    end

    for j = 1:5:N
        plotCellTrajectory!(gaxmain, u, m, j, 3)
    end
    plotCellTrajectory!(gaxmain, u, m, 99, 3)
    #plotCellTrajectory!(gaxmain, u, m,  100, 3)
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotResults1D_spatial_density(u, var)
    #f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    #    resolution=(500, 500))
    f = Figure(fontsize = 32,backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    #gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), xlabel="x", ylabel="y")
    gaxmain = Axis(ga[1, 1], limits=(0, 2*pi, 0, 0.5), xlabel="x", ylabel="y")
    #CRange = findMinMax(var)
    CRange = (0,50)
    for i in eachindex(u)
        if i%5 == 0
            lines!(gaxmain, u[i][:,1], var[i], linewidth=5)
        end
    end
    return f
end

function plotResults1D_Velocity(u, var)
    #f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    #    resolution=(500, 500))
    f = Figure(fontsize = 32,backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    #gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), xlabel="x", ylabel="y")
    gaxmain = Axis(ga[1, 1], limits=(0, 2*pi, 1, 8), aspect=DataAspect(), xlabel="x", ylabel="y")
    #CRange = findMinMax(var)
    CRange = (0,0.4)
    for i in eachindex(u)
        lines!(gaxmain, u[i][:,1], u[i][:,2], color=var[i], colorrange=CRange,
            colormap=:jet, linewidth=5)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=:jet,
        flipaxis=false, label="v [mm/day]") 
    return f
end