using CairoMakie
using ColorSchemes
using Colors

function DiscVSContDensity_plot(Discrete_Solution, Continuum_Solution, index, btype)
    if index > 11
        error("index is too large!")
    end
    # Getting data Continuum
    θ_cont, R_cont, ρ_cont =  Continuum_Solution;
    
    # Getting Data Discrete
    x = zeros(size(Discrete_Solution.u,1),size(Discrete_Solution.u[1],1))
    y = zeros(size(x))
    temp_θ = zeros(size(x))
    temp_ρ = zeros(size(x))
    θ_disc = zeros(size(x))
    ρ_disc = zeros(size(x))
    for i in axes(x,1)
        x[i,:] .= Discrete_Solution.u[i][:,1];
        y[i,:] .= Discrete_Solution.u[i][:,2];
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
        ρ_disc[i,:] = [temp_ρ[i, θmin_idx:end]; temp_ρ[i, 1:θmin_idx-1]]
    end
    # sorting discrete solution values

    
    ## Making Figure
    txtSize = 35;
    tickSize = 25;
    plot_font = "Arial"
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(0, 2π, 22, 60), xticklabelsize = tickSize, yticklabelsize = tickSize, 
                    xlabel="θ", xlabelsize=txtSize, xlabelfont = plot_font,
                    ylabel="ρ", ylabelsize=txtSize,  ylabelfont = plot_font,
                    title = "Pore: $btype, t = $(Discrete_Solution.t[index]) simulation days", titlesize = txtSize, titlefont = plot_font)
    
    # plotting Continuum
    cont_index = 1 + (index - 1)*1000
    cont_line = CairoMakie.lines!(gaxmain, θ_cont, ρ_cont[cont_index,:], linewidth=5, color=:red)


    # plotting Discrete
    disc_index = index;
    disc_scat = CairoMakie.stairs!(gaxmain, θ_disc[disc_index,:], ρ_disc[disc_index,:], step=:pre, linewidth=3, color=:black)
    Legend(f[1,2], [cont_line, disc_scat], ["Continuum", "Discrete m = 1"])

    return f
end