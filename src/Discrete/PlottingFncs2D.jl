
# Colormaps available at: https://docs.juliahub.com/MakieGallery/Ql23q/0.2.17/generated/colors.html#Colormaps

"""
    plotResults2D(u, var, cmap, crange, cbarlabel, D, kf)

Generate a 2D plot to visualize results with lines and scatter points.

# Arguments
- `u::Vector`: A vector of 2D arrays representing the data points.
- `var::Vector`: A vector of values associated with each data point for coloring.
- `cmap::AbstractColorMap`: The colormap used for coloring the plot.
- `crange::AbstractVector`: The color range for mapping values to colors.
- `cbarlabel::AbstractString`: The label for the colorbar.
- `D::Number`: A parameter to be displayed in the plot title.
- `kf::Number`: Another parameter to be displayed in the plot title.

# Returns
- `Figure`: A Makie Figure object representing the 2D plot.
"""
function plotResults2D(u, var, cmap, crange, cbarlabel, D, kf)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), 
              xlabel="x", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="y", ylabelsize = txtSize, yticklabelsize = tickSize,
              title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    for i in eachindex(u)
        plotInterface(gaxmain, u, var, cmap, CRange, i)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotInterface(gaxmain, u, var, cmap, CRange, index)
    CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
            colormap=cmap, linewidth=5)
    CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
        colormap=cmap, markersize=6)
end


function plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, index)
    Lplot = CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
            colormap=cmap, linewidth=5)
    Splot = CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
        colormap=cmap, markersize=6)
    return Lplot,Splot
end


function animateResults2D(u, var, cmap, crange, cbarlabel, D, kf, filename)
    txtSize = 35;
    tickSize = 25;
    CRange = crange
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), 
            xlabel="x", xlabelsize = txtSize, xticklabelsize = tickSize,
            ylabel="y", ylabelsize = txtSize, yticklabelsize = tickSize,
            title = "kₛ = $D, kf = $kf", titlesize = txtSize)
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
            flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, 1)
    Lplot,Splot = plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, 1)

    frames = 2:length(u)+50
    record(f,filename,frames; framerate = 10) do frame
        delete!(f.content[1],Lplot)
        delete!(f.content[1],Splot)
        if frame > length(u)
            index = length(u)
        else
            index = frame
        end
        Lplot,Splot = plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, index)
    end
end



"""
    plotResults2D(u, t, var, cmap, crange, cbarlabel, D, kf)

Generate a 2D plot to visualize results with lines representing angular positions over time.

# Arguments
- `u::Vector`: A vector of 2D arrays representing the data points.
- `t::Vector`: A vector of time values corresponding to the data points.
- `var::Vector`: A vector of values associated with each data point for coloring.
- `cmap::AbstractColorMap`: The colormap used for coloring the plot.
- `crange::AbstractVector`: The color range for mapping values to colors.
- `cbarlabel::AbstractString`: The label for the colorbar.
- `D::Number`: A parameter to be displayed in the plot title.
- `kf::Number`: Another parameter to be displayed in the plot title.

# Returns
- `Figure`: A Makie Figure object representing the 2D plot.
"""
function plotResults2D(u, t, var, cmap, crange, cbarlabel, D, kf)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], 
              xlabel="t", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="θ", ylabelsize = txtSize, yticklabelsize = tickSize,
              title = "kₛ = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    θ = zeros(size(u[1],1)+1,size(t,1))
    ξ = zeros(size(u[1],1)+1,size(t,1))
    for i in eachindex(t)
        x = [u[i][:, 1]; u[i][1,1]].data
        y = [u[i][:, 2]; u[i][1,2]].data
        θ[:,i] = atan.(y,x)
        ξ[:,i] = [var[i]; var[i][1]].data
    end
    for j in axes(θ,1)
        CairoMakie.lines!(gaxmain, t, θ[j,:], color=ξ[j,:], colorrange=CRange,
                colormap=cmap, linewidth=4)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end


## to be fixed
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

