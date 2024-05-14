function CountLengths(Lengths, value)
    return count(value .== Lengths)
end

function plotForceLawCompareStairs(data1)
    function CountLengths(Lengths, value)
        return count(value .== Lengths)
    end
    txtSize = 18;
    tickSize = 18;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(455, 455))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1],
                    xlabel="Length", xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Count", ylabelsize = txtSize, yticklabelsize = tickSize)

    for ii in axes(data1,1)
        data_length = round.(1 ./ data1[ii].data, digits=0)
        lengths = Float64[]
        data1_count = Float64[]

        for length in round.(LinRange(0,21,21), digits=0)
            push!(lengths, length)
            push!(data1_count, CountLengths(data_length, length))
        end
        
        #CairoMakie.barplot!(lengths, data2_count, strokecolor = :black, strokewidth = 1, alpha=0.1)
        #CairoMakie.barplot!(lengths, data1_count, strokecolor = :black, strokewidth = 1, alpha=0.1)
        CairoMakie.stairs!(gaxmain,lengths,data1_count,linewidth=2)
        #CairoMakie.lines!(gaxmain,lengths,data1_count,linewidth=2)

        CairoMakie.xlims!(4, 21)
        CairoMakie.ylims!(0,110)
    end
    return f
end

