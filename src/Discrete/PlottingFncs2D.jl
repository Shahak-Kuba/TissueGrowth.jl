
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
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
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

function plotResults2D(u, var, cmap, crange, cbarlabel, D, kf, axisLims)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-axisLims[1], axisLims[1], -axisLims[2], axisLims[2]), aspect=DataAspect(), 
              xlabel="x [μm]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="y [μm]", ylabelsize = txtSize, yticklabelsize = tickSize)
              #title = "D = $D, kf = $kf", titlesize = txtSize)
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

function plotInterface!(gaxmain, u, var, cmap, CRange, index, lw)
    CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
            colormap=cmap, linewidth=lw)
    CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
        colormap=cmap, markersize=lw+1)
end


function plotResults2D(u, var, cmap, crange, cbarlabel, D, kf, axisLims, embedded_cells, multiInterfaces)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-axisLims[1], axisLims[1], -axisLims[2], axisLims[2]), aspect=DataAspect(), 
              xlabel="x", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="y", ylabelsize = txtSize, yticklabelsize = tickSize,
              title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    if multiInterfaces
        for i in 1:9:size(u,1)
            plotInterface!(gaxmain, u, var, cmap, CRange, i)
        end
    else
        plotInterface!(gaxmain, u, var, cmap, CRange, 1)
        plotInterface!(gaxmain, u, var, cmap, CRange, size(u,1))
    end
    plotEmbeddedCells!(gaxmain, embedded_cells)
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotEmbeddedCells!(gaxmain, embedded_cell_pos)
    for i in axes(embedded_cell_pos,1)
        cell = embedded_cell_pos[i]
        CairoMakie.lines!(gaxmain, cell[1,:], cell[2,:],color=:black,linewidth=8)
    end
end

## AREA COMPARE PLOTTING CODE

# δt compare code

function plotδtAreaResults(Ω₁,t₁,Ω₂,t₂,Ω₃,t₃,N,kf)
    COMPARE = true
    txtSize = 16;
    tickSize = 16;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(455, 400))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], 
              xlabel="t [Days]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="Ω-error [μm^2]", ylabelsize = txtSize, yticklabelsize = tickSize)

    if !COMPARE
        t = LinRange(0,t₁[end],500)
        Ωₐ = Ω_analytic(Ω₁[1],N,kf,t)
        Line0 = plotAreaVsTime!(gaxmain, t, Ωₐ, :green, :solid, "Analytic")
        Line1 = plotAreaVsTime!(gaxmain, t₁, Ω₁, :blue, :solid, "δt = 0.01")
        Line2 = plotAreaVsTime!(gaxmain, t₂, Ω₂, :red, :dash, "δt = 0.001")
        Line3 = plotAreaVsTime!(gaxmain, t₃, Ω₃, :black, :dot, "δt = 0.0001")
        axislegend(gaxmain, merge = true, unique = true)
        #Legend(f[1,2],[Line0,Line1,Line2,Line3], ["Analytic","δt = 0.01", "δt = 0.001","δt = 0.0001"])
    else
        Line1 = plotAreaDiffVsTime!(gaxmain, t₁, Ω₁, N, kf, :blue, :solid, "δt = 0.01")
        Line2 = plotAreaDiffVsTime!(gaxmain, t₂, Ω₂, N, kf, :red, :dash, "δt = 0.001")
        Line3 = plotAreaDiffVsTime!(gaxmain, t₃, Ω₃, N, kf, :black, :dot, "δt = 0.0001")
        #Legend(f[1,2],[Line1,Line2,Line3], ["δt = 0.01", "δt = 0.001","δt = 0.0001"])
        axislegend(gaxmain, merge = true, unique = true, position = :lt)
    end

    return f
end

function plotAreaDiffVsTime!(gaxmain, t, Ωₛ, N, kf, clr, style, name)
    Ωₐ = Ω_analytic(Ωₛ[1],N,kf,t)
    CairoMakie.lines!(gaxmain, t, Ωₛ.-Ωₐ, color=clr, label=name, linewidth=4, linestyle=style)
end


# shape compare plotting code
function plotMultiSimResults2D(Solution, axislims, cmap, CRange)
    txtSize = 16;
    tickSize = 16;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(655, 400))
    ga = f[1, 1] = GridLayout()

    for Diffusivity = axes(Solution,1)
        for Shape = axes(Solution[1],1)
            # Setting gaxmain (axis ticks and labels)
            if Diffusivity == 1
                if Shape == size(Solution[1],1)
                    gaxmain = Axis(ga[Shape, Diffusivity], limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), xticklabelsize = tickSize, yticklabelsize = tickSize)
                else
                    gaxmain = Axis(ga[Shape, Diffusivity], limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize)
                end
            elseif Shape == size(Solution[1],1)
                gaxmain = Axis(ga[Shape, Diffusivity], limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), yticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize)
            else
                gaxmain = Axis(ga[Shape, Diffusivity], limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsvisible = false, yticklabelsize = tickSize)
            end
            # Plotting Interface
            u = Solution[Diffusivity][Shape].u
            var = Solution[Diffusivity][Shape].Density
            for i in eachindex(u)
                plotInterface!(gaxmain, u, var, cmap, CRange, i, 2)
            end
        end
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=15,
        flipaxis=false, label="Density ρ", labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotMultiAreaVsTime(Ω₁,t₁,Ω₂,t₂,N,kf)
    txtSize = 18;
    tickSize = 18;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(455, 455))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], 
                    xlabel="t [Days]", xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Ω [μm^2]", ylabelsize = txtSize, yticklabelsize = tickSize)
    
    t = LinRange(0,t₁[end],500)
    Ωₐ = Ω_analytic(Ω₁[1],N,kf,t)

    Analytic_Sol = plotAreaVsTime!(gaxmain, t, Ωₐ, :red, :solid, "Analytic")
    Square_Sol = plotAreaVsTime!(gaxmain, t₁, Ω₁, :blue, :dash, "Square Pore")
    Hex_Sol = plotAreaVsTime!(gaxmain, t₂, Ω₂, :black, :dot, "Hex Pore")

    #Legend(f[1,1],[Analytic_Sol,Square_Sol,Hex_Sol], ["Analytic Circle", "Discrete Square","Discrete Hex"])
    axislegend(gaxmain, merge = true, unique = true)
    return f
end

function plotAreaVsTime!(gaxmain, t, Ωₛ, clr, style, name)
    CairoMakie.lines!(gaxmain, t, Ωₛ, color=clr, label = name, linewidth=4, linestyle=style)
end

"""
    plotThetaVsTime(u, t, var, cmap, crange, cbarlabel, D, kf)

Generate a plot to visualize results with lines representing angular positions over time.

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
function plotThetaVsTime(u, t, var, cmap, crange, cbarlabel, D, kf)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], 
              xlabel="t [days]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="θ [radians]", ylabelsize = txtSize, yticklabelsize = tickSize)
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



function plotOtValueVsTime(t, Ω, embedded_cell_count, Ot)
    # Sorting Data
    filled_Area = Ω[1] .- Ω
    y = embedded_cell_count./filled_Area
    y[1] = 0.0
    # Creating Figure
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1],
              xlabel="t [days]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="Simulation Ot", ylabelsize = txtSize, yticklabelsize = tickSize,
              title = "Ot = $Ot", titlesize = txtSize)
    
    Ot_line = CairoMakie.lines!(gaxmain, t, Ot.*ones(size(t)), linewidth=3, linestyle = :dash, color = :black)
    Sim_Ot_Line = CairoMakie.lines!(gaxmain, t, y, linewidth=5, color = :blue)
    Legend(f[1,2],[Ot_line,Sim_Ot_Line], ["Ot value", "Simulated Ot"])
    return f
end


# Plot to compare with Buenzli et al. 2020

function plotCompareRegressionBuenzli(Ω_estimate, t, Ωnorm_Analytic, t_Analytic, Ωnorm_Discrete, t_Discrete)
    txtSize = 18;
    tickSize = 18;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(455, 455))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], 
                    xlabel="t [Days]", xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Ω(t)/Ω₀", ylabelsize = txtSize, yticklabelsize = tickSize)
    
    plotAreaVsTime!(gaxmain, t_Discrete, Ωnorm_Discrete, :blue, :solid, "Discrete")
    plotAreaVsTime!(gaxmain, t_Analytic, Ωnorm_Analytic, :red, :dash, "Analytic")
    plotAreaVsTime!(gaxmain, t, Ω_estimate, :black, :dash, "Regression Model")

    axislegend(gaxmain, merge = true, unique = true)
    return f
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