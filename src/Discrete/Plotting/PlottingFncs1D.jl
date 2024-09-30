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
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=650, height=650, limits=(0, 1500, -0.1, 1100), aspect=DataAspect(),
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
    txtSize = 40;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=650, height=650, limits=(0, 1500, -10, 1100), 
              xlabel=L"x \; \text{[μm]}", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel=L"y \; \text{[μm]}", ylabelsize = txtSize, yticklabelsize = tickSize)
    CRange = crange
    for i in eachindex(u)
        plotInterface1D!(gaxmain, u, var, cmap, CRange, i, 7)
    end

    plotCellTrajectory!(gaxmain, u, m, Int(N/4) - 10, 3)
    plotCellTrajectory!(gaxmain, u, m, Int(N/4), 3)
    plotCellTrajectory!(gaxmain, u, m, Int(N/4) + 10, 3)
    #plotCellTrajectory!(gaxmain, u, m,  100, 3)
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,height=650,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

# NoI: Number of Interfaces to plot (equally spaced)

function plotResults1D(u, var, cmap, crange, cbarlabel, D, kf, m, N, NoI)
    txtSize = 45;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=650, height=650, limits=(0, 1500, -10, 2600), 
              xlabel=L"\text{x [μm]}", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel=L"\text{y [μm]}", ylabelsize = txtSize, yticklabelsize = tickSize)
    CRange = crange
    Interface_Step = Int(floor(size(u,1)/NoI))
    for i in 1:Interface_Step:size(u,1)
        plotInterface1D!(gaxmain, u, var, cmap, CRange, i, 7)
    end

    plotCellTrajectory!(gaxmain, u, m, Int(N/4) - 10, 3)
    plotCellTrajectory!(gaxmain, u, m, Int(N/4), 3)
    plotCellTrajectory!(gaxmain, u, m, Int(N/4) + 10, 3)
    plotCellTrajectory!(gaxmain, u, m, Int(N/2), 3)
    plotCellTrajectory!(gaxmain, u, m, 3*Int(N/4) - 10, 3)
    plotCellTrajectory!(gaxmain, u, m, 3*Int(N/4), 3)
    plotCellTrajectory!(gaxmain, u, m, 3*Int(N/4) + 10, 3)
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

function plotThetaVsTime1D(u, t, var, cmap, crange, cbarlabel)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], 
              xlabel="x [mm]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="t [days]", ylabelsize = txtSize, yticklabelsize = tickSize)
    CRange = crange
    x = zeros(size(u[1],1)-1,size(t,1))
    ξ = zeros(size(u[1],1)-1,size(t,1))
    for i in eachindex(t)
        x[:,i] = u[i][2:end, 1].data
        ξ[:,i] = var[i][2:end].data
    end
    for j in axes(x,1)
        CairoMakie.lines!(gaxmain, x[j,:], t, color=ξ[j,:], colorrange=CRange,
                colormap=cmap, linewidth=4)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

# Only middle growth
function plotStationaryBoundary(u, var, cmap, crange, cbarlabel)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-0.01, 1.51, -0.1, 1.1), aspect=DataAspect(), 
              xlabel="x", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="y", ylabelsize = txtSize, yticklabelsize = tickSize)
    CRange = crange
    for index in eachindex(u)
        if typeof(var) == Vector{Vector{Float64}}
            if index == 1
                CairoMakie.lines!(gaxmain, u[index][:,1].data, u[index][:,2].data, color=var[index], colorrange=CRange,
                    colormap=cmap, linewidth=5)
                CairoMakie.scatter!(gaxmain, u[index][:,1].data, u[index][:,2].data, color=var[index], colorrange=CRange,
                    colormap=cmap, markersize=6)
            else
                u_fixed = fixYvalue1D(u[1], u[index])
                CairoMakie.lines!(gaxmain, u_fixed[:,1], u_fixed[:,2], color=var[index], colorrange=CRange,
                    colormap=cmap, linewidth=5)
                CairoMakie.scatter!(gaxmain, u_fixed[:,1], u_fixed[:,2], color=var[index], colorrange=CRange,
                    colormap=cmap, markersize=6)
            end
        else
            if index == 1
                CairoMakie.lines!(gaxmain, u[index][:,1].data, u[index][:,2].data, color=var[index].data, colorrange=CRange,
                    colormap=cmap, linewidth=5)
                CairoMakie.scatter!(gaxmain, u[index][:,1].data, u[index][:,2].data, color=var[index].data, colorrange=CRange,
                    colormap=cmap, markersize=6)
            else
                u_fixed = fixYvalue1D(u[1], u[index])
                CairoMakie.lines!(gaxmain, u_fixed[:,1], u_fixed[:,2], color=var[index].data, colorrange=CRange,
                    colormap=cmap, linewidth=5)
                CairoMakie.scatter!(gaxmain, u_fixed[:,1], u_fixed[:,2], color=var[index].data, colorrange=CRange,
                    colormap=cmap, markersize=6)
            end
        end
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function fixYvalue1D(u0, u1)
    fixed_u1 = zeros(size(u1))
    delta = u1[1,2] - u0[1,2]
    for ii in axes(u1,1)
        fixed_u1[ii,:] .= u1[ii,:] - [0,delta]
    end

    return fixed_u1
end