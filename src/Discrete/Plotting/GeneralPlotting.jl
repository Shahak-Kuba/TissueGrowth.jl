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
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], 
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
    CairoMakie.lines!(gaxmain, hcat(spring_traj...)'[:,1], hcat(spring_traj...)'[:,2], color=:green, linewidth=lw)

end