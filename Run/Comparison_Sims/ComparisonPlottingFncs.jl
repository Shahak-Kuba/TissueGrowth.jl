
function DiscVSContDensity_plot(Discrete_Solution, Continuum_Solution, index)
    if index > 11
        error("index is too large!")
    end
    # Getting data Continuum
    θ_cont, R_cont, ρ_cont =  Continuum_Solution;
    
    # Getting Data Discrete
    x = zeros(size(Discrete_Solution.u,1),size(Discrete_Solution.u[1],1))
    y = zeros(size(x))
    θ_disc = zeros(size(x))
    ρ_disc = zeros(size(x))
    for i in axes(x,1)
        x[i,:] .= Discrete_Solution.u[i][:,1];
        y[i,:] .= Discrete_Solution.u[i][:,2];
        θ_disc[i,:] .= atan.(y[i,:],x[i,:]);
        for j in axes(θ_disc,2)
            if θ_disc[i,j] < 0
                θ_disc[i,j] = θ_disc[i,j] + 2*π
            end
        end
        ρ_disc[i,:] = Discrete_Solution.Density[i].data
    end
    
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
                    title = "t = $(Discrete_Solution.t[index]) simulation days", titlesize = txtSize, titlefont = plot_font)
    
    # plotting Continuum
    cont_index = 1 + (index - 1)*1000
    CairoMakie.lines!(gaxmain, θ_cont, ρ_cont[cont_index,:], linewidth=5, color=:red)


    # plotting Discrete
    disc_index = index;
    CairoMakie.scatter!(gaxmain, θ_disc[disc_index,:], ρ_disc[disc_index,:], markersize=9, color=:black)

    return f
end