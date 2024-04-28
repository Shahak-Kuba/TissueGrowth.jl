
# Colormaps available at: https://docs.juliahub.com/MakieGallery/Ql23q/0.2.17/generated/colors.html#Colormaps


function plotResults2D(u, var, cmap, crange, cbarlabel, axisLims, N, m)
    txtSize = 40;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-axisLims[1], axisLims[1], -axisLims[2], axisLims[2]), aspect=DataAspect(), 
              xlabel=L"\text{x [μm]}", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel=L"\text{y [μm]}", ylabelsize = txtSize, yticklabelsize = tickSize)
              #title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    for i in eachindex(u)
        plotInterface!(gaxmain, u, var, cmap, CRange, i)
    end

    plot_cell_traj = false # User set

    if plot_cell_traj
        for j = 1:3:N
            plotCellTrajectory!(gaxmain, u, m, j, 3)
        end
    end
    #plotCellTrajectory!(gaxmain, u, m, 35, 3)
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotResults2D_Quadrant(u, var, cmap, crange, cbarlabel, axisLims, N, m)
    txtSize = 45;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(0, axisLims[1], 0, axisLims[2]), aspect=DataAspect(), 
              xlabel=L"\text{x [μm]}", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel=L"\text{y [μm]}", ylabelsize = txtSize, yticklabelsize = tickSize)
              #title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    for i in eachindex(u)
        plotInterface!(gaxmain, u, var, cmap, CRange, i, 7)
    end

    for i in 5:5:95
        plotCellTrajectory!(gaxmain, u, m, i, 5)
    end

    #plotCellTrajectory!(gaxmain, u, m, 15, 5)
    #plotCellTrajectory!(gaxmain, u, m, 20, 5)
    #plotCellTrajectory!(gaxmain, u, m, 25, 5)
    #plotCellTrajectory!(gaxmain, u, m, 30, 5)
    #plotCellTrajectory!(gaxmain, u, m, 35, 5)

    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotStress2D_Quadrant(u, var, cmap, Crange, cbarlabel, axisLims)
    txtSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(0, axisLims[1], 0, axisLims[2]), aspect=DataAspect(), 
              xlabel=L"x \; \text{[μm]}", xlabelsize = txtSize+10, xticklabelsize = txtSize,
              ylabel=L"y \; \text{[μm]}", ylabelsize = txtSize+10, yticklabelsize = txtSize)
              #title = "D = $D, kf = $kf", titlesize = txtSize)
    lw = 5
    for index in eachindex(u)
        if index == 1 || index == size(u,1)
            CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=:black,linewidth=lw)
            CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=:black, markersize=lw+1)
        else
            #CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=:grey,linewidth=lw)
            #CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=:grey, markersize=lw+1)    
        end
    end

    for i in axes(u[1],1)
        plotSpringBoundaryTrajectory!(gaxmain, u, var, 5, cmap, Crange, i)
    end

    Colorbar(f[1, 2], limits=Crange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize + 10, ticklabelsize = txtSize)
    return f
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
                    gaxmain = Axis(ga[Shape, Diffusivity], limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), xticks = [-1, 0, 1], xticklabelsize = tickSize, yticklabelsize = tickSize, yticks = [-1, 0, 1])
                else
                    gaxmain = Axis(ga[Shape, Diffusivity], limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), xticks = [-1, 0, 1], xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize, yticks = [-1, 0, 1])
                end
            elseif Shape == size(Solution[1],1)
                gaxmain = Axis(ga[Shape, Diffusivity], limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), yticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize, xticks = [-1, 0, 1], yticks = [-1, 0, 1])
            else
                gaxmain = Axis(ga[Shape, Diffusivity], limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsvisible = false, yticklabelsize = tickSize, xticks = [-1, 0, 1], yticks = [-1, 0, 1])
            end
            # Plotting Interface
            u = Solution[Diffusivity][Shape].u
            var = Solution[Diffusivity][Shape].Vₙ
            #var = Solution[Diffusivity][Shape].Density
            for i in eachindex(u)
                plotInterface!(gaxmain, u, var, cmap, CRange, i, 2)
            end
        end
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=15,
        flipaxis=false, label="Velocity", labelsize = txtSize, ticklabelsize = tickSize)
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

