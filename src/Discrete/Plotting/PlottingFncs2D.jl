
# Colormaps available at: https://docs.juliahub.com/MakieGallery/Ql23q/0.2.17/generated/colors.html#Colormaps


function plotResults2D(u, var, cmap, crange, cbarlabel, axisLims, N, m, NoI)
    txtSize = 35;
    tickSize = 40;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=650, height=650,limits=(-axisLims[1], axisLims[1], -axisLims[2], axisLims[2]), aspect=DataAspect(), 
              xlabel="x [μm]", xlabelsize = txtSize, xticklabelsize = txtSize,
              ylabel="y [μm]", ylabelsize = txtSize, yticklabelsize = txtSize)
              #title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    Interface_Step = Int(floor(size(u,1)/NoI))
    for i in 1:Interface_Step:size(u,1)
        plotInterface!(gaxmain, u, var, cmap, CRange, i, 7)
    end

    plot_cell_traj = false # User set
    if plot_cell_traj
        for j = 1:3:N
            plotCellTrajectory!(gaxmain, u, m, j, 3)
        end
    end
    show_initial_boundaries = false
    if show_initial_boundaries
        #plotting spring boundaries
        CairoMakie.scatter!(gaxmain, [u[1][:, 1]; u[1][1,1]].data, [u[1][:, 2]; u[1][1,2]].data, color="grey", markersize=15)
        #plotting cell boundaries
        for ii in 1:m:size(u[1],1)
            CairoMakie.scatter!(gaxmain, u[1][ii, 1], u[1][ii, 2], color="black", marker=:xcross,markersize=25)
        end
    end

    show_final_boundaries = false
    if show_final_boundaries
        #plotting spring boundaries
        #CairoMakie.scatter!(gaxmain, [u[end][:, 1]; u[end][1,1]].data, [u[end][:, 2]; u[end][1,2]].data, color="grey", markersize=15)
        #plotting cell boundaries
        for ii in 1:m:size(u[1],1)
            CairoMakie.scatter!(gaxmain, u[end][ii, 1], u[end][ii, 2], color="black", marker=:xcross,markersize=25)
        end
    end

    #plotCellTrajectory!(gaxmain, u, m, 35, 3)
    CairoMakie.Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = tickSize, ticklabelsize = txtSize)
    return f
end

function plotResults2D(u, var, cmap, crange, cbarlabel, axisLims, axisTicks, N, m, NoI)
    txtSize = 35;
    tickSize = 40;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=650, height=650, limits=(-axisLims[1], axisLims[1], -axisLims[2], axisLims[2]), aspect=DataAspect(), 
              xlabel="x [μm]", xlabelsize = txtSize, xticklabelsize = tickSize, xticks = axisTicks,
              ylabel="y [μm]", ylabelsize = txtSize, yticklabelsize = tickSize, yticks = axisTicks)
              #title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    Interface_Step = Int(floor(size(u,1)/NoI))
    for i in 1:Interface_Step:size(u,1)
        plotInterface!(gaxmain, u, var, cmap, CRange, i, 6)
    end

    plot_cell_traj = false # User set
    if plot_cell_traj
        for j = 1:3:N
            plotCellTrajectory!(gaxmain, u, m, j, 3)
        end
    end
    show_initial_boundaries = false
    if show_initial_boundaries
        #plotting spring boundaries
        CairoMakie.scatter!(gaxmain, [u[1][:, 1]; u[1][1,1]].data, [u[1][:, 2]; u[1][1,2]].data, color="grey", markersize=15)
        #plotting cell boundaries
        for ii in 1:m:size(u[1],1)
            CairoMakie.scatter!(gaxmain, u[1][ii, 1], u[1][ii, 2], color="black", marker=:xcross,markersize=25)
        end
    end

    show_final_boundaries = false
    if show_final_boundaries
        #plotting spring boundaries
        CairoMakie.scatter!(gaxmain, [u[end][:, 1]; u[end][1,1]].data, [u[end][:, 2]; u[end][1,2]].data, color="grey", markersize=15)
        #plotting cell boundaries
        for ii in 1:m:size(u[1],1)
            CairoMakie.scatter!(gaxmain, u[end][ii, 1], u[end][ii, 2], color="black", marker=:xcross,markersize=25)
        end
    end

    #plotCellTrajectory!(gaxmain, u, m, 35, 3)
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotResults2D(u, var, cmap, crange, cbarlabel, axisLims, xTicks, yTicks, N, m, NoI, 
    show_cell_traj, show_initial_boundaries, show_final_boundaries)
    txtSize = 40;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=650, height=650,limits=(-axisLims[1], axisLims[1], -axisLims[2], axisLims[2]), aspect=DataAspect(), 
              xlabel=L"x\text{ [\mu m]}", xlabelsize = txtSize, xticklabelsize = tickSize, xticks=xTicks,
              ylabel=L"y\text{ [\mu m]}", ylabelsize = txtSize, yticklabelsize = tickSize, yticks=yTicks)
    CRange = crange
    Interface_Step = Int(floor(size(u,1)/NoI))
    for i in 1:Interface_Step:size(u,1)
        plotInterface!(gaxmain, u, var, cmap, CRange, i, 7)
    end

    if show_cell_traj
        for j = 1:3:N
            plotCellTrajectory!(gaxmain, u, m, j, 3)
        end
    end

    if show_initial_boundaries
        #plotting spring boundaries
        CairoMakie.scatter!(gaxmain, [u[1][:, 1]; u[1][1,1]].data, [u[1][:, 2]; u[1][1,2]].data, color="grey", markersize=15)
        #plotting cell boundaries
        for ii in 1:m:size(u[1],1)
            CairoMakie.scatter!(gaxmain, u[1][ii, 1], u[1][ii, 2], color="black", marker=:xcross,markersize=25)
        end
    end

    if show_final_boundaries
        #plotting spring boundaries
        CairoMakie.scatter!(gaxmain, [u[end][:, 1]; u[end][1,1]].data, [u[end][:, 2]; u[end][1,2]].data, color="grey", markersize=15)
        #plotting cell boundaries
        for ii in 1:m:size(u[1],1)
            CairoMakie.scatter!(gaxmain, u[end][ii, 1], u[end][ii, 2], color="black", marker=:xcross,markersize=25)
        end
    end
    
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotResults2D_Quadrant(u, var, cmap, crange, cbarlabel, axisLims, N, m, NoI)
    txtSize = 45;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=650, height=650, limits=(0, axisLims[1], 0, axisLims[2]), aspect=DataAspect(), 
              xlabel="x [μm]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="y [μm]", ylabelsize = txtSize, yticklabelsize = tickSize)
              #title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    Interface_Step = Int(floor(size(u,1)/NoI))
    for i in 1:Interface_Step:size(u,1)
        plotInterface!(gaxmain, u, var, cmap, CRange, i, 7)
    end

    if u[1][1,2] == 0
        for i in 2:5:50
            plotCellTrajectory!(gaxmain, u, m, i, 5)
        end
    else
        for i in 1:5:50
            plotCellTrajectory!(gaxmain, u, m, i, 5)
        end
    end

    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotStress2D_Quadrant(u, var, cmap, Crange, cbarlabel, axisLims)
    txtSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=650, height=650,limits=(0, axisLims[1], 0, axisLims[2]), aspect=DataAspect(), 
              xlabel="x [μm]", xlabelsize = txtSize+10, xticklabelsize = txtSize,
              ylabel="y [μm]", ylabelsize = txtSize+10, yticklabelsize = txtSize)
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


function plotResults2D_embedded(u, var, cmap, crange, cbarlabel, D, kf, axisLims, embedded_cells, multiInterfaces)
    txtSize = 40;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=650, height=650,limits=(-axisLims[1], axisLims[1], -axisLims[2], axisLims[2]), aspect=DataAspect(), 
              xlabel="x [μm]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="y [μm]", ylabelsize = txtSize, yticklabelsize = tickSize)
              #title = "D = $D, kf = $kf", titlesize = txtSize)
    CRange = crange
    if multiInterfaces
        for i in 1:20:size(u,1)
            plotInterface!(gaxmain, u, var, cmap, CRange, i)
        end
        plotInterface!(gaxmain, u, var, cmap, CRange, size(u,1))
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
    txtSize = 24;
    tickSize = 24;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(650, 600))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], 
              xlabel="t [Days]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="Rel-Error (Ω)", ylabelsize = txtSize, yticklabelsize = tickSize)

    if !COMPARE
        t = LinRange(0,t₁[end],500)
        Ωₐ = Ω_analytic(Ω₁[1],N,kf,t)
        Line0 = plotAreaVsTime!(gaxmain, t, Ωₐ, :green, :solid, "Analytic")
        Line1 = plotAreaVsTime!(gaxmain, t₁, Ω₁, :blue, :solid, "δt = 0.01")
        Line2 = plotAreaVsTime!(gaxmain, t₂, Ω₂, :red, :dash, "δt = 0.001")
        Line3 = plotAreaVsTime!(gaxmain, t₃, Ω₃, :black, :dot, "δt = 0.0001")
        #axislegend(gaxmain, merge = true, unique = true)
        #Legend(f[1,2],[Line0,Line1,Line2,Line3], ["Analytic","δt = 0.01", "δt = 0.001","δt = 0.0001"])
    else
        Line1 = plotAreaRelErrorVsTime!(gaxmain, t₁, Ω₁, N, kf, :blue, :solid, "Δt = 1")
        Line2 = plotAreaRelErrorVsTime!(gaxmain, t₂, Ω₂, N, kf, :red, :solid, "Δt = 0.01")
        Line3 = plotAreaRelErrorVsTime!(gaxmain, t₃, Ω₃, N, kf, :black, :solid, "Δt = 0.0001")
        #Legend(f[1,2],[Line1,Line2,Line3], ["δt = 0.01", "δt = 0.001","δt = 0.0001"])
        axislegend(gaxmain, merge = true, unique = true, position = :lt)
    end

    return f
end

function plotAreaDiffVsTime!(gaxmain, t, Ωₛ, N, kf, clr, style, name)
    Ωₐ = Ω_analytic(Ωₛ[1],N,kf,t)
    CairoMakie.lines!(gaxmain, t, Ωₐ.-Ωₛ, color=clr, label=name, linewidth=4, linestyle=style)
end

function plotAreaRelErrorVsTime!(gaxmain, t, Ωₛ, N, kf, clr, style, name)
    Ωₐ = Ω_analytic(Ωₛ[1],N,kf,t)
    CairoMakie.lines!(gaxmain, t, ((Ωₛ.-Ωₐ)./Ωₐ), color=clr, label=name, linewidth=4, linestyle=style)
end


# shape compare plotting code
function plotMultiSimResults2D(Solution, axislims, cmap, CRange)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()

    for Diffusivity = axes(Solution,1)
        for Shape = axes(Solution[1],1)
            # Setting gaxmain (axis ticks and labels)
            if Diffusivity == 1
                if Shape == size(Solution[1],1)
                    gaxmain = Axis(ga[Shape, Diffusivity], height = 650, width=650, limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), xticks = [-1, 0, 1], xticklabelsize = tickSize, yticklabelsize = tickSize, yticks = [-1, 0, 1])
                else
                    gaxmain = Axis(ga[Shape, Diffusivity], height = 650, width=650, limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), xticks = [-1, 0, 1], xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize, yticks = [-1, 0, 1])
                end
            elseif Shape == size(Solution[1],1)
                gaxmain = Axis(ga[Shape, Diffusivity], height = 650, width=650, limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), yticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize, xticks = [-1, 0, 1], yticks = [-1, 0, 1])
            else
                gaxmain = Axis(ga[Shape, Diffusivity], height = 650, width=650, limits=(-axislims[1], axislims[1], -axislims[2], axislims[2]), xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsvisible = false, yticklabelsize = tickSize, xticks = [-1, 0, 1], yticks = [-1, 0, 1])
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

function plotMultiAreaVsTime(t_discrete_large, t_discrete_mid, t_discrete_small, Ω_large_square, Ω_large_hex, Ω_mid_square, Ω_mid_hex, Ω_small_square, Ω_small_hex, N_large, N_mid, N_small, kf)
    txtSize = 40;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(900, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 650, width=650, limits=(0,48,0,1250000),
                    xlabel="t [Days]", xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Ω [μm²]", ylabelsize = txtSize, yticklabelsize = tickSize)
    
    # Small
    t = LinRange(0,t_discrete_small[end],500)
    Ωₐ = Ω_analytic(Ω_small_square[1],N_small,kf,t)
    Analytic_Sol = plotAreaVsTime!(gaxmain, t, Ωₐ, :red, :solid, L"\text{Analytic}", 7)
    Square_Sol_ = plotAreaVsTime!(gaxmain, t_discrete_small, Ω_small_square, :black, :dash, L"\text{Square Pore}", 6)
    Hex_Sol_Hook = plotAreaVsTime!(gaxmain, t_discrete_small, Ω_small_hex, :blue, :dot, L"\text{Hex Pore}", 6)

    # Medium
    t = LinRange(0,t_discrete_mid[end],500)
    Ωₐ = Ω_analytic(Ω_mid_square[1],N_mid,kf,t)
    Analytic_Sol = plotAreaVsTime!(gaxmain, t, Ωₐ, :red, :solid, 7)
    Square_Sol_ = plotAreaVsTime!(gaxmain, t_discrete_mid, Ω_mid_square, :black, :dash, 6)
    Hex_Sol_Hook = plotAreaVsTime!(gaxmain, t_discrete_mid, Ω_mid_hex, :blue, :dot, 6)

    # Large
    t = LinRange(0,t_discrete_large[end],500)
    Ωₐ = Ω_analytic(Ω_large_square[1],N_large,kf,t)
    Analytic_Sol = plotAreaVsTime!(gaxmain, t, Ωₐ, :red, :solid, 7)
    Square_Sol_ = plotAreaVsTime!(gaxmain, t_discrete_large, Ω_large_square, :black, :dash, 6)
    Hex_Sol_Hook = plotAreaVsTime!(gaxmain, t_discrete_large, Ω_large_hex, :blue, :dot, 6)


    #Legend(f[1,1],[Analytic_Sol,Square_Sol,Hex_Sol], ["Analytic Circle", "Discrete Square","Discrete Hex"])
    axislegend(gaxmain, merge = true, unique = true, labelsize=txtSize)
    return f
end

function plotMultiAreaVsTime(t_discrete, Ω_large_square, Ω_large_hex, Ω_mid_square, Ω_mid_hex, Ω_small_square, Ω_small_hex, N_large, N_mid, N_small, kf, Ω_hex_ρ₀_small, Ω_hex_ρ₀_mid, Ω_hex_ρ₀_large)
    txtSize = 40;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(900, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 650, width=650, limits=(0,24,0,1400000),
                    xlabel="t [Days]", xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Ω [μm²]", ylabelsize = txtSize, yticklabelsize = tickSize, yticks=[250000, 560000, 1000000])
    
    t = LinRange(0,t_discrete[end],500)

    # Small
    Ωₐ = Ω_analytic(Ω_small_square[1],N_small,kf,t)
    Analytic_Sol = plotAreaVsTime!(gaxmain, t, Ωₐ, :red, :solid, L"\text{Analytic}", 7)
    Square_Sol_ = plotAreaVsTime!(gaxmain, t_discrete, Ω_small_square, :black, :dash, L"\text{Square Pore}", 6)
    Hex_Sol_Hook = plotAreaVsTime!(gaxmain, t_discrete, Ω_small_hex, :blue, :dot, L"\text{Hex Pore}", 6)
    Hex_Sol_ρ₀ = plotAreaVsTime!(gaxmain, t_discrete, Ω_hex_ρ₀_small, :green, :solid, L"\text{Hex Pore q_{0}=0.05}", 6)


    # Medium
    Ωₐ = Ω_analytic(Ω_mid_square[1],N_mid,kf,t)
    Analytic_Sol = plotAreaVsTime!(gaxmain, t, Ωₐ, :red, :solid, 7)
    Square_Sol_ = plotAreaVsTime!(gaxmain, t_discrete, Ω_mid_square, :black, :dash, 6)
    Hex_Sol_Hook = plotAreaVsTime!(gaxmain, t_discrete, Ω_mid_hex, :blue, :dot, 6)
    Hex_Sol_ρ₀ = plotAreaVsTime!(gaxmain, t_discrete, Ω_hex_ρ₀_mid, :green, :solid, 6)


    # Large
    Ωₐ = Ω_analytic(Ω_large_square[1],N_large,kf,t)
    Analytic_Sol = plotAreaVsTime!(gaxmain, t, Ωₐ, :red, :solid, 7)
    Square_Sol_ = plotAreaVsTime!(gaxmain, t_discrete, Ω_large_square, :black, :dash, 6)
    Hex_Sol_Hook = plotAreaVsTime!(gaxmain, t_discrete, Ω_large_hex, :blue, :dot, 6)
    Hex_Sol_ρ₀ = plotAreaVsTime!(gaxmain, t_discrete, Ω_hex_ρ₀_large, :green, :solid, 6)


    #Legend(f[1,1],[Analytic_Sol,Square_Sol,Hex_Sol], ["Analytic Circle", "Discrete Square","Discrete Hex"])
    axislegend(gaxmain, merge = true, unique = true, labelsize=txtSize)
    return f
end

#function plotAreaVsTime!(gaxmain, t, Ωₛ, clr, style, name)
#    CairoMakie.lines!(gaxmain, t, Ωₛ, color=clr, label = name, linewidth=4, linestyle=style)
#end
function plotAreaVsTime!(gaxmain, t, Ωₛ, clr, style, lw)
    CairoMakie.lines!(gaxmain, t, Ωₛ, color=clr, linewidth=lw, linestyle=style)
end

function plotAreaVsTime!(gaxmain, t, Ωₛ, clr, style, name, lw)
    CairoMakie.lines!(gaxmain, t, Ωₛ, color=clr, label = name, linewidth=lw, linestyle=style)
end



# Plot to compare with Buenzli et al. 2020

function plotCompareRegressionBuenzli(Ω_estimate, t, Ωnorm_Analytic, t_Analytic, Ωnorm_Discrete, t_Discrete)
    txtSize = 40;
    tickSize = 35;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(850, 850))
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

