function plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, index)
    Lplot = CairoMakie.lines!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
            colormap=cmap, linewidth=8)
    Splot = CairoMakie.scatter!(gaxmain, [u[index][:, 1]; u[index][1,1]].data, [u[index][:, 2]; u[index][1,2]].data, color=[var[index]; var[index][1]].data, colorrange=CRange,
        colormap=cmap, markersize=9)
    return Lplot,Splot
end


function animateResults2D(t, u, var, cmap, crange, cbarlabel, filename)
    txtSize = 45;
    tickSize = 40;
    CRange = crange
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-200, 200, -200, 200), aspect=DataAspect(), 
            xlabel="x [μm]", xlabelsize = txtSize, xticklabelsize = tickSize,
            ylabel="y [μm]", ylabelsize = txtSize, yticklabelsize = tickSize,
            title = "t = $(t[1])days", titlesize = txtSize)
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
            flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, 1)
    Lplot,Splot = plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, 1)

    frames = 2:length(u)+50
    record(f,filename,frames; framerate = 10) do frame
        delete!(f.content[1],Lplot)
        delete!(f.content[1],Splot)
        if frame > length(u)
            index = length(u)
        else
            index = frame
        end
        T = round(t[index];digits=2)
        gaxmain.title="t = $T days"
        Lplot,Splot = plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, index)
    end
end

function animateResults2D(t, u, var, cmap, crange, cbarlabel, filename, embedded_cell_pos, embed_times)
    txtSize = 45;
    tickSize = 40;
    CRange = crange
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(-200, 200, -200, 200), aspect=DataAspect(), 
            xlabel="x [μm]", xlabelsize = txtSize, xticklabelsize = tickSize,
            ylabel="y [μm]", ylabelsize = txtSize, yticklabelsize = tickSize,
            title = "t = $(t[1])days", titlesize = txtSize)
    Colorbar(f[1, 2], limits=CRange, colormap=cmap, size=30,
            flipaxis=false, label=cbarlabel, labelsize = txtSize, ticklabelsize = tickSize)
    plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, 1)
    Lplot,Splot = plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, 1)
    reverse!(embedded_cell_pos)

    frames = 2:length(u)+50
    record(f,filename,frames; framerate = 10) do frame
        delete!(f.content[1],Lplot)
        delete!(f.content[1],Splot)
        if frame > length(u)
            index = length(u)
        else
            index = frame
        end
        T = round(t[index];digits=2)
        gaxmain.title="t = $T days"
        Lplot,Splot = plotInterfaceAnimation(gaxmain, u, var, cmap, CRange, index)
        for ii in eachindex(embed_times)
            if t[index] >= embed_times[ii]
                CairoMakie.lines!(gaxmain, embedded_cell_pos[ii][1,:], embedded_cell_pos[ii][2,:],color=:black,linewidth=8)
            end
        end
    end
end