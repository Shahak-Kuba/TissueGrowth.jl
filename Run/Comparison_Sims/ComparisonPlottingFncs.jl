using CairoMakie
using ColorSchemes
using Colors

function DiscVSContDensity_plot(gaxmain, Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, index)
    if index > 11
        error("index is too large!")
    end
    # Getting data Continuum
    θ_cont, R_cont, ρ_cont =  Continuum_Solution;
    max_y = ceil(maximum(ρ_cont[end,:]); sigdigits=1)
    if max_y - maximum(ρ_cont[end,:]) > 5
        max_y = ceil(maximum(ρ_cont[end,:]); sigdigits=1)+5
    end
    # Getting Data Discrete
    θ_disc1, R_disc1, ρ_disc1 = Convert_Discrete_Data(Discrete_Solution_m1, m1)
    θ_disc2, R_disc2, ρ_disc2 = Convert_Discrete_Data(Discrete_Solution_m2, m2)
   

    # plotting Discrete
    disc_index = index;
    disc_stair1 = CairoMakie.stairs!(gaxmain, θ_disc1[disc_index,:], ρ_disc1[disc_index,:], step=:center, linewidth=3, color=:blue)
    disc_stair2 = CairoMakie.stairs!(gaxmain, θ_disc2[disc_index,:], ρ_disc2[disc_index,:], step=:center, linewidth=3, color=:green)

    # plotting Continuum
    cont_index = 1 + (index - 1)*1000
    cont_line = CairoMakie.lines!(gaxmain, θ_cont, ρ_cont[cont_index,:], linewidth=2, color=:red)

    text!(gaxmain, 0.2, max_y-1.2 ,text= "t=$(Discrete_Solution_m1.t[index])", fontsize=24)
    
    return cont_line, disc_stair1, disc_stair2
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

function DiscVSContDensity_plot_all(Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies, num_cols)
     # Getting data Continuum
    θ_cont, R_cont, ρ_cont =  Continuum_Solution;
    min_y = floor(ρ_cont[1,1]; sigdigits=1)
    max_y = ceil(maximum(ρ_cont[end,:]); sigdigits=1)
    max_y = ceil(maximum(ρ_cont[end,:]); sigdigits=1)
    if max_y - maximum(ρ_cont[end,:]) > 5
        max_y = ceil(maximum(ρ_cont[end,:]); sigdigits=1)+5
    end
    
    txtSize = 35;
    tickSize = 25;
    plot_font = "Arial"
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1200, 1000))
    ga = f[1, 1] = GridLayout()
    col_size = num_cols;
    row_size = Int64(ceil(length(indicies)/col_size))
    for i in axes(indicies,1)
        row = Int64(floor((i-1)/col_size) + 1)
        col = Int64((i-1)%col_size + 1)
        # only showing needed axis ticks
        if col == 1
            if row == row_size
                gaxmain = Axis(ga[row, col], limits=(0, 2π, min_y, max_y), xticklabelsize = tickSize, yticklabelsize = tickSize)
            else
                gaxmain = Axis(ga[row, col], limits=(0, 2π, min_y, max_y), xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize)
            end
        elseif row == row_size
            gaxmain = Axis(ga[row, col], limits=(0, 2π, min_y, max_y), yticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize)
        else
            gaxmain = Axis(ga[row, col], limits=(0, 2π, min_y, max_y), xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsvisible = false, yticklabelsize = tickSize)
        end
        cont_line, disc_stair1, disc_stair2 = DiscVSContDensity_plot(gaxmain, Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies[i])
        Legend(f[1,2], [cont_line, disc_stair1, disc_stair2], ["Continuum", "Discrete m = $m1", "Discrete m = $m2"], labelsize=tickSize)
    end
    Label(ga[0, :], "Pore: Square", fontsize = 45)
    Label(ga[:, 0], "Density ρ", fontsize = 30, rotation=π/2)
    Label(ga[row_size+1, :], "Angular position θ", fontsize = 30)
    return f
end

function DiscVSContShape_plot(Discrete_Solution, m, Continuum_Solution, xbound, ybound, cmap, Cbar_min, Cbar_max)
    # getting Continuum data
    θ_cont, R_cont, ρ_cont =  Continuum_Solution;
    # getting Discrete data
    x_disc = zeros(size(Discrete_Solution.u,1),size(Discrete_Solution.u[1],1)+1)
    y_disc = zeros(size(x_disc))
    ρ_disc = zeros(size(x_disc))
    for i in axes(x_disc,1)
        x_disc[i,:] .= [Discrete_Solution.u[i][:,1];Discrete_Solution.u[i][1,1]];
        y_disc[i,:] .= [Discrete_Solution.u[i][:,2];Discrete_Solution.u[i][1,2]];
        ρ_disc[i,:] .= [Discrete_Solution.Density[i].data;Discrete_Solution.Density[i].data[1]]./m;
    end

    txtSize = 35;
    tickSize = 20;
    plot_font = "Arial"
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(600, 800))
    ga = f[1, 1] = GridLayout()
    col_size = 1;
    #row_size = Int64(ceil(length(indicies)/col_size))
    Cbar_range = (Cbar_min, Cbar_max)
    
    # plotting Discrete
    gaxmain = Axis(ga[1, 1], limits=(-xbound, xbound, -ybound, ybound), xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize)
    for i in axes(x_disc,1)
        CairoMakie.lines!(gaxmain, x_disc[i,:],  y_disc[i,:], color=ρ_disc[i,:], colorrange=Cbar_range, colormap=cmap, linewidth=5)
        CairoMakie.scatter!(gaxmain, x_disc[i,:],  y_disc[i,:], color=ρ_disc[i,:], colorrange=Cbar_range, colormap=cmap, markersize=6)
    end
    # plotting Continuum
    gbxmain = Axis(ga[2, 1], limits=(-xbound, xbound, -ybound, ybound), xticklabelsize = tickSize, yticklabelsize = tickSize)
    for i in 1:1000:size(R_cont,1)
        CairoMakie.lines!(gbxmain, [R_cont[i,:]; R_cont[i,1]].*cos.([θ_cont;θ_cont[1]]), [R_cont[i,:]; R_cont[i,1]].*sin.([θ_cont;θ_cont[1]]), color=[ρ_cont[i,:];ρ_cont[i,1]], colorrange=Cbar_range,
            colormap=cmap, linewidth=5)
        CairoMakie.scatter!(gbxmain, [R_cont[i,:]; R_cont[i,1]].*cos.([θ_cont;θ_cont[1]]), [R_cont[i,:]; R_cont[i,1]].*sin.([θ_cont;θ_cont[1]]), color=[ρ_cont[i,:];ρ_cont[i,1]], colorrange=Cbar_range,
            colormap=cmap, markersize=6)
    end
    Colorbar(f[1, 2], limits=Cbar_range, size=20, ticklabelsize = tickSize, colormap=cmap,
        flipaxis=false, label="Density ρ [cells/length]", labelsize=txtSize)

    return f
end