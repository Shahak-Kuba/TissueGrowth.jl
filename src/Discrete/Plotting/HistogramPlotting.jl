function plotAttributeHistogram(data,Label::String)
    txtSize = 18;
    tickSize = 18;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(455, 455))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], 
                    xlabel=Label, xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Count", ylabelsize = txtSize, yticklabelsize = tickSize)
    clr = :red
    plotAttributeHistogram!(gaxmain, data, clr)
    return f
end

function plotAttributeHistogram!(gaxmain, data, clr)
    hist_data = round.(data,digits=4)
    Makie.hist!(gaxmain,hist_data, normalize=:density,bins = 3,strokewidth = 1, strokecolor = :black)
end

function plotForceLawCompareHistogram(data1, data2, Label)
    txtSize = 18;
    tickSize = 18;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(455, 455))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1],
                    xlabel=Label, xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Count", ylabelsize = txtSize, yticklabelsize = tickSize)
    clr1 = :red
    plotAttributeHistogram!(gaxmain, data1, clr1)
    clr2 = :blue
    plotAttributeHistogram!(gaxmain, data2, clr2)
    CairoMakie.xlims!(5, 20)
    CairoMakie.ylims!(0,110)
    return f
end