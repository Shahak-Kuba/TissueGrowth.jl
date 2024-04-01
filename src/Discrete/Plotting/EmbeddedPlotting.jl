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

function plotEmbeddedCells!(gaxmain, embedded_cell_pos)
    for i in axes(embedded_cell_pos,1)
        cell = embedded_cell_pos[i]
        CairoMakie.lines!(gaxmain, cell[1,:], cell[2,:],color=:black,linewidth=8)
    end
end