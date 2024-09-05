"""
    plotThetaVsTime(u, t, var, cmap, crange, cbarlabel, D, kf)

Generate a plot to visualize results with lines representing angular positions over time.

# Arguments
- `u::Vector`: A vector of 2D arrays representing the data points.
- `t::Vector`: A vector of time values corresponding to the data points.
- `var::Vector`: A vector of values associated with each data point for coloring.
- `cmap::AbstractColorMap`: The colormap used for coloring the plot.
- `crange::AbstractVector`: The color range for mapping values to colors.
- `cbarlabel::AbstractString`: The label for the colorbar.
- `D::Number`: A parameter to be displayed in the plot title.
- `kf::Number`: Another parameter to be displayed in the plot title.

# Returns
- `Figure`: A Makie Figure object representing the 2D plot.
"""
function plotThetaVsTime(u, t, var, cmap, crange, cbarlabel)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1, 1, 1),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 650, width=650,
              xlabel="t [days]", xlabelsize = txtSize, xticklabelsize = tickSize,
              ylabel="θ [radians]", ylabelsize = txtSize, yticklabelsize = tickSize)
    CRange = crange
    θ = zeros(size(u[1],1)+1,size(t,1))
    ξ = zeros(size(u[1],1)+1,size(t,1))
    for i in eachindex(t)
        x = [u[i][:, 1]; u[i][1,1]].data
        y = [u[i][:, 2]; u[i][1,2]].data
        θ[:,i] = atan.(y,x)
        ξ[:,i] = [var[i]; var[i][1]].data
    end
    for j in axes(θ,1)
        CairoMakie.lines!(gaxmain, t, θ[j,:], color=ξ[j,:], colorrange=CRange,
            colormap=cmap, linewidth=4)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotThetaVsTime_Quadrant(u, t, var, cmap, crange, cbarlabel)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1, 1, 1),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 650, width=650,
              xlabel="θ [radians]", xlabelsize = txtSize, xticklabelsize = tickSize, xticks = ([0, π/4, π/2],["0", "π/4", "π/2"]),
              ylabel="t [days]", ylabelsize = txtSize, yticklabelsize = tickSize)
    CRange = crange
    θ = zeros(size(u[1],1),size(t,1))
    ξ = zeros(size(u[1],1),size(t,1))
    for i in eachindex(t)
        x = u[i][:, 1].data
        y = u[i][:, 2].data
        θ[:,i] = atan.(y,x)
        ξ[:,i] = var[i].data
    end
    QuadSize = Int(size(θ,1)/4)
    for j in Int(QuadSize/2)+1:Int(3*QuadSize/2)+1
        CairoMakie.lines!(gaxmain, θ[j,:], t, color=ξ[j,:], colorrange=CRange,
            colormap=cmap, linewidth=4)
    end
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size = 30,
        flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    return f
end

function plotInterface!(gaxmain, u, var, cmap, CRange, index)
    CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
            colormap=cmap, linewidth=5)
    CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
        colormap=cmap, markersize=6)
end

function plotInterface!(gaxmain, u, var, cmap, CRange, index, lw)
    if typeof(var) == Vector{Vector{Float64}}
        CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]], colorrange=CRange,
            colormap=cmap, linewidth=lw)
        CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]], colorrange=CRange,
            colormap=cmap, markersize=lw+1)
    else
        CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
                colormap=cmap, linewidth=lw)
        CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
            colormap=cmap, markersize=lw+1)
    end
end

function plotInterface_Quadrant!(gaxmain, u, var, cmap, CRange, index, lw)
    start_index = Int(size(u[index],1)/8)
    end_index = Int(3*size(u[index],1)/8)+1

    if typeof(var) == Vector{Vector{Float64}}
        CairoMakie.lines!(gaxmain, u[index][start_index:end_index, 1].data, u[index][start_index:end_index, 2].data, color=var[index][start_index:end_index], colorrange=CRange,
            colormap=cmap, linewidth=lw)
        CairoMakie.scatter!(gaxmain, u[index][start_index:end_index, 1].data, u[index][start_index:end_index, 2].data, color=var[index][start_index:end_index], colorrange=CRange,
            colormap=cmap, markersize=lw+1)
    else
        CairoMakie.lines!(gaxmain, u[index][start_index:end_index, 1].data, u[index][start_index:end_index, 2].data, color=var[index][start_index:end_index].data, colorrange=CRange,
                colormap=cmap, linewidth=lw)
        CairoMakie.scatter!(gaxmain, u[index][start_index:end_index, 1].data, u[index][start_index:end_index, 2].data, color=var[index][start_index:end_index].data, colorrange=CRange,
            colormap=cmap, markersize=lw+1)
    end
end

function plotInterface1D!(gaxmain, u, var, cmap, CRange, index, lw)
    if typeof(var) == Vector{Vector{Float64}}
        CairoMakie.lines!(gaxmain, u[index][:, 1].data, u[index][:, 2].data, color=var[index], colorrange=CRange,
            colormap=cmap, linewidth=lw)
        CairoMakie.scatter!(gaxmain, u[index][:, 1].data, u[index][:, 2].data, color=var[index], colorrange=CRange,
            colormap=cmap, markersize=lw+1)
    else
        CairoMakie.lines!(gaxmain, u[index][:, 1].data, u[index][:, 2].data, color=var[index].data, colorrange=CRange,
                colormap=cmap, linewidth=lw)
        CairoMakie.scatter!(gaxmain, u[index][:, 1].data, u[index][:, 2].data, color=var[index].data, colorrange=CRange,
            colormap=cmap, markersize=lw+1)
    end
end

function plotCellTrajectory!(gaxmain, u, m, cell_index, lw)
    left_cell_boundary_idx = cell_index*m - (m-1)
    right_cell_boundary_idx = left_cell_boundary_idx + m

    left_cell_traj = []
    spring_traj = []
    right_cell_traj = []

    for ii in axes(u,1)
        push!(left_cell_traj, u[ii][left_cell_boundary_idx,:])
        push!(right_cell_traj, u[ii][right_cell_boundary_idx,:])
        push!(spring_traj, u[ii][left_cell_boundary_idx+1:right_cell_boundary_idx-1,:]')
    end

    CairoMakie.lines!(gaxmain, hcat(left_cell_traj...)'[:,1], hcat(left_cell_traj...)'[:,2], color=:black, linewidth=lw)
    #CairoMakie.arrows(gaxmain, hcat(left_cell_traj...)'[end-1,1], hcat(left_cell_traj...)'[end-1,2], hcat(left_cell_traj...)'[end,1], hcat(left_cell_traj...)'[end,2], color=:black, arrowsize=10)
    CairoMakie.lines!(gaxmain, hcat(right_cell_traj...)'[:,1], hcat(right_cell_traj...)'[:,2], color=:black, linewidth=lw)
    CairoMakie.lines!(gaxmain, hcat(spring_traj...)'[:,1], hcat(spring_traj...)'[:,2], color=:red, linewidth=lw)

end

function plotSpringBoundaryTrajectory!(gaxmain, u, var, lw, cmap, Crange, idx)

        spring_boundary_traj = []
        spring_boundary_var = []
        
        for ii in axes(u,1)
            push!(spring_boundary_traj, u[ii][idx,:])
            push!(spring_boundary_var, var[ii][idx])
        end

        CairoMakie.lines!(gaxmain, hcat(spring_boundary_traj...)'[:,1], hcat(spring_boundary_traj...)'[:,2], color=vcat(spring_boundary_var...), colorrange=Crange,colormap=cmap, linewidth=lw)

end

function plotForceLawCompareStairs(data1, t, binSize)
    function CountLengths(Lengths, value)
        return count(value .== Lengths)
    end
    txtSize = 40;
    tickSize = 40;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 900))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], width=850, height=850,
                    xlabel="Length", xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Count", ylabelsize = txtSize, yticklabelsize = tickSize)
    colors = [:darkorange, :red, :green, :purple, :blue]
    for ii in axes(data1,1)
        data_length = round.(1 ./ data1[ii].data, digits=1)
        lengths = Float64[]
        data1_count = Float64[]

        for length in 4:0.1:22
            push!(lengths, length)
            push!(data1_count, CountLengths(data_length, length))
        end

        # new array based on user specified bin size
        data_count_bin_unit = Float64[]
        data_count_bin_3 = Float64[]
        lengths2 = Float64[]

        for jj in 1:10:171
            count_of_unit_length = sum(data1_count[jj:jj+9])
            push!(data_count_bin_unit,count_of_unit_length)
            push!(lengths2, lengths[jj])
        end

        lengths = Float64[3,3,5,5,8,8,11,11,14,14,17,17,20,20,23,23]
        push!(data_count_bin_3, 0)
        push!(data_count_bin_3, 0)
        for jj in 1:3:16
            push!(data_count_bin_3, sum(data_count_bin_unit[jj:jj+2]))
            push!(data_count_bin_3, sum(data_count_bin_unit[jj:jj+2]))
        end
        push!(data_count_bin_3, 0)
        push!(data_count_bin_3, 0)

        time = t[ii]
        
        #CairoMakie.barplot!(lengths, data2_count, strokecolor = :black, strokewidth = 1, alpha=0.1)
        #CairoMakie.barplot!(lengths, data1_count, strokecolor = :black, strokewidth = 1, alpha=0.1)
        if binSize == 3
            CairoMakie.stairs!(gaxmain,lengths,data_count_bin_3,linewidth=4, color= colors[ii], label="t = $time")
        else
            CairoMakie.stairs!(gaxmain,lengths2.+0.5,data_count_bin_unit,linewidth=4, color= colors[ii], label="t = $time")
        end
        #CairoMakie.lines!(gaxmain,lengths,data1_count,linewidth=2)

        CairoMakie.xlims!(0, 25)
        CairoMakie.ylims!(0,110)
    end
    axislegend(gaxmain, merge = true, unique = true, labelsize=tickSize, position = :lt)
    return f
end