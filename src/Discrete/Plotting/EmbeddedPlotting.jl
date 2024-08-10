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
              xlabel="t [days]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="Ot [#/μm²]", ylabelsize = txtSize, yticklabelsize = tickSize)
              #title = "Ot = $Ot", titlesize = txtSize)
    
    Ot_line = CairoMakie.lines!(gaxmain, t, Ot.*ones(size(t)), linewidth=3, linestyle = :dash, color = :black, label = "Expected")
    Sim_Ot_Line = CairoMakie.lines!(gaxmain, t, y, linewidth=5, color = :red, label = "Simulated")
    #Legend(f[1,1],[Ot_line,Sim_Ot_Line], ["Ot value", "Simulated Ot"])
    CairoMakie.xlims!(gaxmain,(0,t[end]))
    axislegend(gaxmain, merge = true, unique = true, labelsize=txtSize)
    return f
end

function plotOtValueVsTime(t, numerical_Ot, set_Ot, min_numerical_Ot, max_numerical_Ot, m)
    # Creating Figure
    txtSize = 45;
    tickSize = 40;
    f = Figure(backgroundcolor=RGBf(1, 1, 1),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 650, width=650,
              xlabel="t [days]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="Ot [#/μm²]", ylabelsize = txtSize, yticklabelsize = tickSize)
              #title = "Ot = $set_Ot", titlesize = txtSize)
    
    Sim_Range = CairoMakie.band!(gaxmain, t, min_numerical_Ot, max_numerical_Ot, color=(:blue,0.2))
    Ot_line = CairoMakie.lines!(gaxmain, t, set_Ot.*ones(size(t)), linewidth=3, linestyle = :dash, color = :black, label = "Ot value")
    Sim_Ot_Line = CairoMakie.lines!(gaxmain, t, numerical_Ot, linewidth=5, color = :red, label = "Simulated Ot")
    CairoMakie.xlims!(gaxmain,(0,t[end]))
    CairoMakie.ylims!(gaxmain,(0,2*set_Ot))
    #Legend(f[1,2],[Ot_line,Sim_Ot_Line], ["Ot value", "Simulated Ot"])
    axislegend(gaxmain, merge = true, unique = true, labelsize=txtSize)
    return f
end

function plotOtValueVsTime(t, numerical_Ot, set_Ot, min_numerical_Ot, max_numerical_Ot, m, std_Ot)
    # Creating Figure
    txtSize = 45;
    tickSize = 40;
    f = Figure(backgroundcolor=RGBf(1, 1, 1),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 600, width=600,
              xlabel="t [days]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="Ot [#/μm²]", ylabelsize = txtSize, yticklabelsize = tickSize)
              #title = "Ot = $set_Ot", titlesize = txtSize)
    
    #Sim_Range = CairoMakie.band!(gaxmain, t, min_numerical_Ot, max_numerical_Ot, color=(:blue,0.2))
    Ot_line = CairoMakie.lines!(gaxmain, t, set_Ot.*ones(size(t)), linewidth=3, linestyle = :dash, color = :black, label = "Expected")
    Sim_Ot_Line = CairoMakie.lines!(gaxmain, t, numerical_Ot, linewidth=5, color = :red, label = "Simulated")
    std_upper = CairoMakie.lines!(gaxmain, t, numerical_Ot .+ std_Ot, linewidth=5, color = :black, label = "Std dev")
    std_lower = CairoMakie.lines!(gaxmain, t, numerical_Ot .- std_Ot, linewidth=5, color = :black)
    CairoMakie.xlims!(gaxmain,(0,t[end]))
    CairoMakie.ylims!(gaxmain,(0,2*set_Ot))
    #Legend(f[1,2],[Ot_line,Sim_Ot_Line], ["Ot value", "Simulated Ot"])
    axislegend(gaxmain, merge = true, unique = true, labelsize=35)
    return f
end

function plotEmbeddedCells!(gaxmain, embedded_cell_pos)
    for i in axes(embedded_cell_pos,1)
        cell = embedded_cell_pos[i]
        CairoMakie.lines!(gaxmain, cell[1,:], cell[2,:],color=:black,linewidth=8)
    end
end

function animateOt(t, Ot, exptected_Ot, filename)
    txtSize = 45;
    tickSize = 40;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 1000))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(0, t[end], 0, 2*exptected_Ot), 
            xlabel="t [days]", xlabelsize = txtSize, xticklabelsize = tickSize,
            ylabel="Ot [#/μm²]", ylabelsize = txtSize, yticklabelsize = tickSize,
            title = "t = $(t[1])days", titlesize = txtSize)

    CairoMakie.lines!(gaxmain, t,exptected_Ot.*ones(size(t)),color=:black, linewidth=7, linestyle = :dash)

    frames = 2:length(t)+50
    record(f,filename,frames; framerate = 10) do frame
        if frame > length(t)
            index = length(t)
        else
            index = frame
        end
        T = round(t[index];digits=2)
        gaxmain.title="t = $T days"
        CairoMakie.lines!(gaxmain, t[1:index],Ot[1:index],color=:red, linewidth=7)
    end
end