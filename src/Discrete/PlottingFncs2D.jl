
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
        plotInterface!(gaxmain, u, var, cmap, CRange, i)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotInterface!(gaxmain, u, var, cmap, CRange, index)
    CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
            colormap=cmap, linewidth=5)
    CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
        colormap=cmap, markersize=6)
end

function plotEmbeddedCells!()
    # TODO
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



function plotEmbeddedCellsVsTime(Embedded_Cells, m, t, Ot)
    # Sorting Data
    
    # Creating Figure
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), 
              xlabel="x", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="y", ylabelsize = txtSize, yticklabelsize = tickSize,
              title = "D = $D, kf = $kf", titlesize = txtSize)
end


# ANIMATION CODE

function plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, index)
    Lplot = CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
            colormap=cmap, linewidth=5)
    Splot = CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
        colormap=cmap, markersize=6)
    return Lplot,Splot
end


function animateResults2D(t, u, var, cmap, crange, cbarlabel, D, kf, filename)
    txtSize = 35;
    tickSize = 25;
    CRange = crange
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-1.5, 1.5, -1.5, 1.5), aspect=DataAspect(), 
            xlabel="x", xlabelsize = txtSize, xticklabelsize = tickSize,
            ylabel="y", ylabelsize = txtSize, yticklabelsize = tickSize,
            title = "t = $(t[1])", titlesize = txtSize)
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
        T = round(t[index];digits=2)
        gaxmain.title="t = $T"
        Lplot,Splot = plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, index)
    end
end