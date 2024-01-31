using CairoMakie
using ColorSchemes
using Colors

function DiscVSContDensity_plot(Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, index, btype)
    if index > 11
        error("index is too large!")
    end
    # Getting data Continuum
    θ_cont, R_cont, ρ_cont =  Continuum_Solution;
    
    # Getting Data Discrete
    θ_disc1, R_disc1, ρ_disc1 = Convert_Discrete_Data(Discrete_Solution_m1, m1)
    θ_disc2, R_disc2, ρ_disc2 = Convert_Discrete_Data(Discrete_Solution_m2, m2)
    # sorting discrete solution values

    min_y = floor(ρ_cont[1,1]; sigdigits=1)
    max_y = ceil(maximum(ρ_cont[end,:]); sigdigits=1)

    
    ## Making Figure
    txtSize = 35;
    tickSize = 25;
    plot_font = "Arial"
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(0, 2π, min_y, max_y), xticklabelsize = tickSize, yticklabelsize = tickSize, 
                    xlabel="θ", xlabelsize=txtSize, xlabelfont = plot_font,
                    ylabel="ρ", ylabelsize=txtSize,  ylabelfont = plot_font,
                    title = "Pore: $btype, t = $(Discrete_Solution.t[index]) Simulation days", titlesize = txtSize, titlefont = plot_font)
    
    # plotting Continuum
    cont_index = 1 + (index - 1)*1000
    cont_line = CairoMakie.lines!(gaxmain, θ_cont, ρ_cont[cont_index,:], linewidth=5, color=:red)


    # plotting Discrete
    disc_index = index;
    disc_scat1 = CairoMakie.stairs!(gaxmain, θ_disc1[disc_index,:], ρ_disc1[disc_index,:], step=:pre, linewidth=3, color=:blue)
    disc_scat2 = CairoMakie.stairs!(gaxmain, θ_disc2[disc_index,:], ρ_disc2[disc_index,:], step=:pre, linewidth=3, color=:black)
    #Legend(f[1,2], [cont_line, disc_scat1, disc_scat2], ["Continuum", "Discrete m = $m1", "Discrete m = $m2"])
    axislegend(gaxmain,[cont_line, disc_scat1, disc_scat2], ["Continuum", "Discrete m = $m1", "Discrete m = $m2"],position=:rt)
    return f
end

function Convert_Discrete_Data(Discrete_Solution, m)
    # Getting Data Discrete
    x = zeros(size(Discrete_Solution.u,1),size(Discrete_Solution.u[1],1))
    y = zeros(size(x))
    temp_r = zeros(size(x))
    temp_θ = zeros(size(x))
    temp_ρ = zeros(size(x))
    θ_disc = zeros(size(x))
    ρ_disc = zeros(size(x))
    r_disc = zeros(size(x))
    for i in axes(x,1)
        x[i,:] .= Discrete_Solution.u[i][:,1];
        y[i,:] .= Discrete_Solution.u[i][:,2];
        temp_r[i,:] .= sqrt.(x[i,:].^2 .+ y[i,:].^2)
        temp_θ[i,:] .= atan.(y[i,:],x[i,:]);
        for j in axes(temp_θ,2)
            if temp_θ[i,j] < 0
                temp_θ[i,j] = temp_θ[i,j] + 2*π
            end
        end
        # reorganising the data from θmin -> θmax
        θmin_idx = argmin(temp_θ[i,:]);
        θ_disc[i,:] = [temp_θ[i, θmin_idx:end]; temp_θ[i, 1:θmin_idx-1]]

        temp_ρ[i,:] = Discrete_Solution.Density[i].data
        ρ_disc[i,:] = [temp_ρ[i, θmin_idx:end]; temp_ρ[i, 1:θmin_idx-1]]./m
        r_disc[i,:] = [temp_r[i, θmin_idx:end]; temp_r[i, 1:θmin_idx-1]]
    end

    return θ_disc,r_disc,ρ_disc
end