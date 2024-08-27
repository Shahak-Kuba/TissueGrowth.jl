using CairoMakie
using ColorSchemes
using Colors

function DiscVSContDensity_plot(gaxmain, Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, index)
    if index > 11
        error("index is too large!")
    end
    # Getting data Continuum
    θ_cont, R_cont, ρ_cont =  Continuum_Solution;
    # Getting Data Discrete
    θ_disc1, R_disc1, ρ_disc1 = Convert_Discrete_Data(Discrete_Solution_m1, m1)
    θ_disc2, R_disc2, ρ_disc2 = Convert_Discrete_Data(Discrete_Solution_m2, m2)
   
    
    # plotting Discrete
    disc_index = index;
    disc_stair1 = CairoMakie.stairs!(gaxmain, θ_disc1[disc_index,:], ρ_disc1[disc_index,:], step=:center, linewidth=7, color=:blue)
    disc_stair2 = CairoMakie.stairs!(gaxmain, θ_disc2[disc_index,:], ρ_disc2[disc_index,:], step=:center, linewidth=7, color=:green)

    # plotting Continuum
    cont_index = 1 + (index - 1)*1000
    cont_line = CairoMakie.lines!(gaxmain, θ_cont, ρ_cont[cont_index,:], linewidth=6, color=:red, linestyle=:solid)
    
    return cont_line, disc_stair1, disc_stair2
end

function DiscVSContDensity_plot(gaxmain, Discrete_Solution_m, m, Continuum_Solution, index)
    if index > 11
        error("index is too large!")
    end
    # Getting data Continuum
    θ_cont, R_cont, ρ_cont =  Continuum_Solution;
    # Getting Data Discrete
    θ_disc, R_disc, ρ_disc = Convert_Discrete_Data(Discrete_Solution_m, m)
   
    if m == 1
        clr = :blue
    else
        clr = :blue
    end
    
    # plotting Continuum
    cont_index = 1 + (index - 1)*1000
    cont_line = CairoMakie.lines!(gaxmain, θ_cont, ρ_cont[cont_index,:], linewidth=3, color=:red, linestyle=:solid)

    # plotting Discrete
    disc_index = index;
    disc_stair = CairoMakie.stairs!(gaxmain, [0 ; θ_disc[disc_index,:];2π], [ρ_cont[cont_index,1]; ρ_disc[disc_index,:]; ρ_cont[cont_index,end]], step=:center, linewidth=3, color=:blue)
    
    return cont_line, disc_stair
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
        ρ_disc[i,:] = [temp_ρ[i, θmin_idx:end]; temp_ρ[i, 1:θmin_idx-1]]
        r_disc[i,:] = [temp_r[i, θmin_idx:end]; temp_r[i, 1:θmin_idx-1]]
    end

    return θ_disc,r_disc,ρ_disc
end

function DiscVSContDensity_plot_all(Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies, num_cols)
    max_y = 0.18
    
    txtSize = 16;
    tickSize = 14;
    plot_font = "Arial"
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(655, 400))
    ga = f[1, 1] = GridLayout()
    col_size = num_cols;
    row_size = Int64(ceil(length(indicies)/col_size))
    for i in axes(indicies,1)
        row = Int64(floor((i-1)/col_size) + 1)
        col = Int64((i-1)%col_size + 1)
        # only showing needed axis ticks
        if col == 1
            clr = :blue
            if row == row_size
                gaxmain = Axis(ga[row, col], limits=(0, 2π, 0, max_y), xticks = ([0, π/2, π, 3π/2, 2π],["0", "π/2", "π", "3π/2", "2π"]), xticklabelsize = tickSize, yticklabelsize = tickSize, yticks = [0, 0.04, 0.08, 0.12, 0.16])
            else
                gaxmain = Axis(ga[row, col], limits=(0, 2π, 0, max_y), xticks = ([0, π/2, π, 3π/2, 2π],["0", "π/2", "π", "3π/2", "2π"]), xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize, yticks = [0, 0.04, 0.08, 0.12, 0.16])
            end
                DiscVSContDensity_plot(gaxmain, Discrete_Solution_m1, m1, Continuum_Solution, indicies[i])
        elseif row == row_size
            clr = :green
            gaxmain = Axis(ga[row, col], limits=(0, 2π, 0, max_y), xticks = ([0, π/2, π, 3π/2, 2π],["0", "π/2", "π", "3π/2", "2π"]),  yticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsize = tickSize, yticks = [0, 0.04, 0.08, 0.12, 0.16])
            DiscVSContDensity_plot(gaxmain, Discrete_Solution_m2, m2, Continuum_Solution, indicies[i])
        else
            clr = :green
            gaxmain = Axis(ga[row, col], limits=(0, 2π, 0, max_y), xticks = ([0, π/2, π, 3π/2, 2π],["0", "π/2", "π", "3π/2", "2π"]), xticklabelsvisible = false, xticklabelsize = tickSize, yticklabelsvisible = false, yticklabelsize = tickSize, yticks = [0, 0.04, 0.08, 0.12, 0.16])
            DiscVSContDensity_plot(gaxmain, Discrete_Solution_m2, m2, Continuum_Solution, indicies[i])
        end
        #cont_line, disc_stair1, disc_stair2 = DiscVSContDensity_plot(gaxmain, Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies[i])
        #Legend(f[1,1], [cont_line, disc_stair1, disc_stair2], ["Continuum", "Discrete m = $m1", "Discrete m = $m2"], labelsize=tickSize)
    end
    #Label(ga[0, :], "Pore: Square", fontsize = 45)
    Label(ga[:, 0], L"\text{Density} \; q \; \text{[1/μm]}", fontsize = 16, rotation=π/2)
    Label(ga[row_size+1, :], L"\text{Angle} \; θ", fontsize = 16)
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
        ρ_disc[i,:] .= [Discrete_Solution.Density[i].data;Discrete_Solution.Density[i].data[1]];
    end

    txtSize = 35;
    tickSize = 28;
    plot_font = "Arial"
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(700, 1000))
    ga = f[1, 1] = GridLayout()
    col_size = 1;
    #row_size = Int64(ceil(length(indicies)/col_size))
    Cbar_range = (Cbar_min, Cbar_max)
    
    # plotting Discrete
    gaxmain = Axis(ga[1, 1], width=450, height=450,limits=(-xbound, xbound, -ybound, ybound), xticklabelsvisible = false, 
     xticklabelsize = tickSize, yticklabelsize = tickSize)
    for i in axes(x_disc,1)
        #if i == size(x_disc,1)
            
        #else
            CairoMakie.lines!(gaxmain, x_disc[i,:],  y_disc[i,:], color=ρ_disc[i,:], colorrange=Cbar_range, colormap=cmap, linewidth=4)
            CairoMakie.scatter!(gaxmain, x_disc[i,:],  y_disc[i,:], color=ρ_disc[i,:], colorrange=Cbar_range, colormap=cmap, markersize=5)
        #end
    end
    
    # plotting Continuum
    gbxmain = Axis(ga[2, 1], width=450, height=450,limits=(-xbound, xbound, -ybound, ybound), xticklabelsize = tickSize, yticklabelsize = tickSize)
    for i in 1:1000:size(R_cont,1)
        CairoMakie.lines!(gbxmain, [R_cont[i,:]; R_cont[i,1]].*cos.([θ_cont;θ_cont[1]]), [R_cont[i,:]; R_cont[i,1]].*sin.([θ_cont;θ_cont[1]]), color=[ρ_cont[i,:];ρ_cont[i,1]], colorrange=Cbar_range,
            colormap=cmap, linewidth=4)
        CairoMakie.scatter!(gbxmain, [R_cont[i,:]; R_cont[i,1]].*cos.([θ_cont;θ_cont[1]]), [R_cont[i,:]; R_cont[i,1]].*sin.([θ_cont;θ_cont[1]]), color=[ρ_cont[i,:];ρ_cont[i,1]], colorrange=Cbar_range,
            colormap=cmap, markersize=5)
    end
    Colorbar(f[1, 2], limits=Cbar_range, size=30, ticklabelsize = tickSize, colormap=cmap,
        flipaxis=false, label=L"\text{Density} \; q \; \text{[1/length]}", labelsize=txtSize)

    return f
end


function DiscVSContShape_plot_all(Discrete_Solutions, m, Continuum_Solutions, xbound, ybound, cmap, Cbar_min, Cbar_max)
    txtSize = 35;
    tickSize = 28;
    plot_font = "Arial"
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(2100, 1000))
    #ga = f[1, 1] = GridLayout()
    ga = f
    Cbar_range = (Cbar_min, Cbar_max)


    for ii in axes(Discrete_Solutions,1)
        # getting Continuum data
        θ_cont, R_cont, ρ_cont =  Continuum_Solutions[ii];
        # getting Discrete data
        x_disc = zeros(size(Discrete_Solutions[ii].u,1),size(Discrete_Solutions[ii].u[1],1)+1)
        y_disc = zeros(size(x_disc))
        ρ_disc = zeros(size(x_disc))
        for i in axes(x_disc,1)
            x_disc[i,:] .= [Discrete_Solutions[ii].u[i][:,1];Discrete_Solutions[ii].u[i][1,1]];
            y_disc[i,:] .= [Discrete_Solutions[ii].u[i][:,2];Discrete_Solutions[ii].u[i][1,2]];
            ρ_disc[i,:] .= [Discrete_Solutions[ii].Density[i].data;Discrete_Solutions[ii].Density[i].data[1]];
        end

        # plotting Discrete
        if ii == 1
            gaxmain = Axis(ga[1, ii], width=450, height=450,limits=(-xbound, xbound, -ybound, ybound), xticks = ([-50, -25, 0, 25, 50]),xticklabelsvisible = false, 
                xticklabelsize = tickSize, yticks = ([-50, -25, 0, 25, 50]), yticklabelsize = tickSize, title="", titlesize = 50)
        else
            gaxmain = Axis(ga[1, ii], width=450, height=450,limits=(-xbound, xbound, -ybound, ybound), xticks = ([-50, -25, 0, 25, 50]),xticklabelsvisible = false, 
            xticklabelsize = tickSize, yticklabelsize = tickSize, yticks = ([-50, -25, 0, 25, 50]), yticklabelsvisible=false, title="", titlesize = 50)
        end
        for i in axes(x_disc,1)
            #if i == size(x_disc,1)
                
            #else
                CairoMakie.lines!(gaxmain, x_disc[i,:],  y_disc[i,:], color=ρ_disc[i,:], colorrange=Cbar_range, colormap=cmap, linewidth=4)
                CairoMakie.scatter!(gaxmain, x_disc[i,:],  y_disc[i,:], color=ρ_disc[i,:], colorrange=Cbar_range, colormap=cmap, markersize=5)
            #end
        end
        
        # plotting Continuum
        if ii == 1
            gbxmain = Axis(ga[2, ii], width=450, height=450,limits=(-xbound, xbound, -ybound, ybound), xticks = ([-50, -25, 0, 25, 50]), xticklabelsize = tickSize, yticks = ([-50, -25, 0, 25, 50]), yticklabelsize = tickSize, title="",titlesize = 50)
        else
            gbxmain = Axis(ga[2, ii], width=450, height=450,limits=(-xbound, xbound, -ybound, ybound), xticks = ([-50, -25, 0, 25, 50]), xticklabelsize = tickSize, yticks = ([-50, -25, 0, 25, 50]), yticklabelsize = tickSize, yticklabelsvisible=false, title="",titlesize = 50)
        end
        for i in 1:1000:size(R_cont,1)
            CairoMakie.lines!(gbxmain, [R_cont[i,:]; R_cont[i,1]].*cos.([θ_cont;θ_cont[1]]), [R_cont[i,:]; R_cont[i,1]].*sin.([θ_cont;θ_cont[1]]), color=[ρ_cont[i,:];ρ_cont[i,1]], colorrange=Cbar_range,
                colormap=cmap, linewidth=4)
            CairoMakie.scatter!(gbxmain, [R_cont[i,:]; R_cont[i,1]].*cos.([θ_cont;θ_cont[1]]), [R_cont[i,:]; R_cont[i,1]].*sin.([θ_cont;θ_cont[1]]), color=[ρ_cont[i,:];ρ_cont[i,1]], colorrange=Cbar_range,
                colormap=cmap, markersize=5)
        end
    end

    Colorbar(f[1:2, 4], limits=Cbar_range, size=30, ticklabelsize = tickSize, colormap=cmap,
            flipaxis=false, label=L"\text{Density} \; q \; \text{[1/length]}", labelsize=txtSize)
    Makie.resize_to_layout!(f)
    return f
end