
using CairoMakie
using ColorSchemes
using Colors
# Colormaps available at: https://docs.juliahub.com/MakieGallery/Ql23q/0.2.17/generated/colors.html#Colormaps

function plotResults2D(u, var, cmap)
    #f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    #    resolution=(500, 500))
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), xlabel="x", ylabel="y")
    #gaxmain = Axis(ga[1, 1], limits=(0, 2*pi, 1, 8), aspect=DataAspect(), xlabel="x", ylabel="y")
    #CRange = findMinMax(var)
    CRange = (0,120)
    for i in eachindex(u)
        lines!(gaxmain, [u[i][:, 1]; u[i][1,1]], [u[i][:, 2]; u[i][1,2]], color=[var[i]; var[i][1]], colorrange=CRange,
            colormap=cmap, linewidth=5)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap,
        flipaxis=false, label="ρ [μm²]")
    return f
end

function plotResults2D_Velocity(u, var, cmap)
    #f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    #    resolution=(500, 500))
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), xlabel="x", ylabel="y")
    #gaxmain = Axis(ga[1, 1], limits=(0, 2*pi, 1, 8), aspect=DataAspect(), xlabel="x", ylabel="y")
    #CRange = findMinMax(var)
    CRange = (0,0.05)
    for i in eachindex(u)
        lines!(gaxmain, [u[i][:, 1]; u[i][1,1]], [u[i][:, 2]; u[i][1,2]], color=[var[i]; var[i][1]], colorrange=CRange,
            colormap=cmap, linewidth=5)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap,
        flipaxis=false, label="Vₙ [μm/s]")
    return f
end


function plotResults1D(u, var)
    #f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    #    resolution=(500, 500))
    f = Figure(fontsize = 32,backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    #gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), xlabel="x", ylabel="y")
    gaxmain = Axis(ga[1, 1], limits=(0, 2*pi, 1, 8), aspect=DataAspect(), xlabel="x", ylabel="y")
    CRange = findMinMax(var)
    #CRange = (0.05,0.2)
    for i in eachindex(u)
        lines!(gaxmain, u[i][:,1], u[i][:,2], color=var[i].data, colorrange=CRange,
            colormap=:jet, linewidth=5)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=:jet,
        flipaxis=false, label="v [mm/day]") 
    return f
end


function plotAreaVStime(sols)
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(700, 500))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], title="Void area over time", xlabel="Time [Days]", ylabel="Ω [μm²]")
    lin1 = lines!(gaxmain, sols[1].t, sols[1].Ω, linewidth=4, linestyle=:solid)
    lin2 = lines!(gaxmain, sols[2].t, sols[2].Ω, linewidth=4, linestyle=:dash)
    lin3 = lines!(gaxmain, sols[3].t, sols[3].Ω, linewidth=4, linestyle=:dot)
    lin4 = lines!(gaxmain, sols[4].t, sols[4].Ω, linewidth=4, linestyle=:dashdot)
    lin5 = lines!(gaxmain, sols[5].t, sols[5].Ω, linewidth=4, linestyle=:dashdot)
    Legend(f[1, 2], [lin1, lin2, lin3, lin4, lin5], ["circle", "triangle", "square", "hex", "star"])
    return f
end

function plotAreaVStime_δt(sol_δt1, sol_δt2, sol_δt3, sol_δt4, Geometry)
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(700, 500))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], title="Void area difference compared with solution for δt = 0.00001 for $Geometry pore", xlabel="Time [Days]", ylabel="Ω-difference [units²]")
    lin1 = lines!(gaxmain, sol_δt2.t, sol_δt2.Ω - sol_δt1.Ω, linewidth=4, linestyle=:solid)
    lin2 = lines!(gaxmain, sol_δt3.t, sol_δt3.Ω - sol_δt1.Ω, linewidth=4, linestyle=:dash)
    lin3 = lines!(gaxmain, sol_δt4.t, sol_δt4.Ω - sol_δt1.Ω, linewidth=4, linestyle=:dot)
    Legend(f[1, 2], [lin1, lin2, lin3], ["δt = 0.0001", "δt = 0.001", "δt = 0.01"])
    return f
end

function plotKapVsVel(sol)
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        resolution=(700, 500))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], title="Curvature vs velocity @ different days", xlabel="κ", ylabel="vₙ [μms⁻¹]")
    for ii in axes(sol.Κ,1)
        if sol.btype == "triangle"
            mid = floor(Int64, length(sol.Κ[ii])/3)
        elseif sol.btype == "square"
            mid = floor(Int64, length(sol.Κ[ii])/4)
        elseif sol.btype == "hex"
            mid = floor(Int64, length(sol.Κ[ii])/6)
        elseif sol.btype == "circle"
            mid = floor(Int64, length(sol.Κ[ii])/2)
        elseif sol.btype == "star"
            mid = floor(Int64, length(sol.Κ[ii])/1)
        end

        x = sort(sol.Κ[ii][2:mid])
        y = sort(sol.Vₙ[ii][2:mid])

        scatter!(gaxmain, x, y, markersize = 25, marker = '*')
        

        #fit = curve_fit(LogFit, x, y)
        #yfit = fit.(x)
        #lines!(gaxmain, x, yfit, linewidth=3, linestyle=:dash)
    end

    return f
end

