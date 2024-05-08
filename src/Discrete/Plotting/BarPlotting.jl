function CountLengths(Lengths, value)
    return count(value .== Lengths)
end

function plotForceLawCompareBarplot(data1, data2)
    lengths = Int64[]
    data1_count = Int64[]
    data2_count = Int64[]

    for length in 5:20
        push!(lengths, length)
        push!(data1_count, CountLengths(data1, length))
        push!(data2_count, CountLengths(data2, length))
    end

    txtSize = 18;
    tickSize = 18;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(455, 455))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1],
                    xlabel="Length", xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Count", ylabelsize = txtSize, yticklabelsize = tickSize)
    
    CairoMakie.barplot!(lengths, data2_count, strokecolor = :black, strokewidth = 1, alpha=0.1)
    CairoMakie.barplot!(lengths, data1_count, strokecolor = :black, strokewidth = 1, alpha=0.1)

    CairoMakie.xlims!(4, 21)
    CairoMakie.ylims!(0,100)
    return f
end