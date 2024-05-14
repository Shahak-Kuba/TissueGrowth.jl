 
function plotContinuumResults_Cartesian(x, h, ρ, cmap, D, kf)
    txtSize = 35;
    tickSize = 25;
    plot_font = "Arial"
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(0, 2π, 1, 8), aspect=DataAspect(), xticklabelsize = tickSize, yticklabelsize = tickSize, 
                    xlabel="x", xlabelsize=txtSize, xlabelfont = plot_font,
                    ylabel="h(x,t)", ylabelsize=txtSize,  ylabelfont = plot_font,
                    title = "D = $D, kf = $kf", titlesize = txtSize, titlefont = plot_font)
    CRange = (20,50)
    for i in 1:150:size(h,1)
        lines!(gaxmain, x, h[i,:], color=ρ[i,:], colorrange=CRange,
            colormap=cmap, linewidth=5)
        scatter!(gaxmain, x, h[i,:], color=ρ[i,:], colorrange=CRange,
            colormap=cmap,markersize = 6)
    end
    Colorbar(f[1, 2], limits=CRange, size=20, ticklabelsize = txtSize, colormap=cmap,
        flipaxis=false, label="Density ρ [cells/length]", labelsize=txtSize)
    return f
end


function plotContinuumResults_Polar(θ, R, ρ, cmap, D, kf)
    txtSize = 35;
    tickSize = 25;
    plot_font = "Arial"
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
        size=(1000, 800))
    ga = f[1, 1] = GridLayout()
    lim = 60
    gaxmain = Axis(ga[1, 1], limits=(-lim, lim, -lim, lim), aspect=DataAspect(), 
                    xlabel="x", xlabelsize=txtSize, xlabelfont = plot_font,
                    ylabel="y", ylabelsize=txtSize,  ylabelfont = plot_font,
                    title = "D = $D, kf = $kf", titlesize = txtSize, titlefont = plot_font)
    CRange = (0.05,0.2)
    for i in 1:1000:size(R,1)
        CairoMakie.lines!(gaxmain, [R[i,:]; R[i,1]].*cos.([θ;θ[1]]), [R[i,:]; R[i,1]].*sin.([θ;θ[1]]), color=[ρ[i,:];ρ[i,1]], colorrange=CRange,
            colormap=cmap, linewidth=5)
        CairoMakie.scatter!(gaxmain, [R[i,:]; R[i,1]].*cos.([θ;θ[1]]), [R[i,:]; R[i,1]].*sin.([θ;θ[1]]), color=[ρ[i,:];ρ[i,1]], colorrange=CRange,
            colormap=cmap, markersize=6)
    end
    Colorbar(f[1, 2], limits=CRange, size=20, ticklabelsize = txtSize, colormap=cmap,
        flipaxis=false, label="Density ρ [cells/length]", labelsize=txtSize)
    return f
end