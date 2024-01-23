using CairoMakie
using ColorSchemes
using Colors

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

function plotResults1D(u,var,D,kf,cmap,upLim,lowLim)
    txtSize = 35;
    tickSize = 25;
    plot_font = "Arial"
    f = Figure(fontsize = 32,backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(0, 2π, 1, 8), aspect=DataAspect(), xticklabelsize = tickSize, yticklabelsize = tickSize, 
                    xlabel="x", xlabelsize=txtSize, xlabelfont = plot_font,
                    ylabel="y", ylabelsize=txtSize,  ylabelfont = plot_font,
                    title = "D = $D, kf = $kf", titlesize = txtSize, titlefont = plot_font)
    CRange = (upLim,lowLim)
    for i in eachindex(u)
        lines!(gaxmain, [u[i][:,1]; u[i][1,1]], [u[i][:,2]; u[i][1,2]], color=[var[i].data; var[i].data[1]], colorrange=CRange,
            colormap=cmap, linewidth=5)
        scatter!(gaxmain, [u[i][:,1]; u[i][1,1]], [u[i][:,2]; u[i][1,2]], color=[var[i].data; var[i].data[1]], colorrange=CRange,
            colormap=cmap,markersize = 6)
        #lines!(gaxmain, u[i][1,:], u[i][2,:], linewidth=5)
    end
    Colorbar(f[1, 2], limits=CRange, size=20, ticklabelsize = txtSize, colormap=cmap,
        flipaxis=false, label="Density ρ [cells/length]", labelsize=txtSize)
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