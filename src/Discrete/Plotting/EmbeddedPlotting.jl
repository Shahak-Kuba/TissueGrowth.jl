function plotOtValueVsTime(t, Ω, embedded_cell_count, Ot, m)
    # Sorting Data
    filled_Area = Ω[1] .- Ω
    y = embedded_cell_count./filled_Area
    y[1] = 0.0
    # Creating Figure
    txtSize = 35;
    tickSize = 30;
    f = Figure(backgroundcolor=RGBf(1, 1, 1),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 650, width=650,
              xlabel=L"t \; \text{[days]}", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel=L"Ot\;[#/\mu\text{m}^{2}]", ylabelsize = txtSize, yticklabelsize = tickSize)
              #title = "Ot = $Ot", titlesize = txtSize)
    
    Ot_line = CairoMakie.lines!(gaxmain, t, Ot.*ones(size(t)).*m, linewidth=3, linestyle = :dash, color = :black, label = L"\text{Ot value}")
    Sim_Ot_Line = CairoMakie.lines!(gaxmain, t, y.*m, linewidth=5, color = :blue, label = L"\text{Simulated Ot}")
    #Legend(f[1,1],[Ot_line,Sim_Ot_Line], ["Ot value", "Simulated Ot"])
    axislegend(gaxmain, merge = true, unique = true, labelsize=txtSize)
    return f
end

function plotOtValueVsTime(t, numerical_Ot, set_Ot, min_numerical_Ot, max_numerical_Ot)
    # Creating Figure
    txtSize = 35;
    tickSize = 30;
    f = Figure(backgroundcolor=RGBf(1, 1, 1),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 650, width=650,
              xlabel="t [days]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="Simulation Ot", ylabelsize = txtSize, yticklabelsize = tickSize,
              title = "Ot = $set_Ot", titlesize = txtSize)
    
    Sim_Range = CairoMakie.band!(gaxmain, t, min_numerical_Ot, max_numerical_Ot, color=(:blue,0.2))
    Ot_line = CairoMakie.lines!(gaxmain, t, set_Ot.*ones(size(t)), linewidth=3, linestyle = :dash, color = :black)
    Sim_Ot_Line = CairoMakie.lines!(gaxmain, t, numerical_Ot, linewidth=5, color = :red)
    CairoMakie.xlims!(gaxmain,(0,t[end]))
    CairoMakie.ylims!(gaxmain,(0,2*set_Ot))
    Legend(f[1,2],[Ot_line,Sim_Ot_Line], ["Ot value", "Simulated Ot"])
    
    return f
end

function plotEmbeddedCells!(gaxmain, embedded_cell_pos)
    for i in axes(embedded_cell_pos,1)
        cell = embedded_cell_pos[i]
        CairoMakie.lines!(gaxmain, cell[1,:], cell[2,:],color=:black,linewidth=8)
    end
end