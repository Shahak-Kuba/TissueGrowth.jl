using TissueGrowth
using Makie
using CairoMakie


## Resting length investigation (with a₁ and a₂ + no phase shift)
Aₕ = (l₀, l₁, p) -> p.k .* ((l₁.^2 - l₀.^2)./2 + p.a .* (l₀ - l₁))
Aₙ = (l₀, l₁, p) -> p.k .* p.a^2 .* ((l₁ .- l₀)./p.a .+ log.(l₀) .- log.(l₁))

l_min = 5
l_max = 20.0
l = LinRange(l_min, l_max, 100)

ks = 7.5
a_nonlinear = 5.0
p_nonlinear = (k = ks, a = a_nonlinear)

a_hookean = ((l_max - l_min)^2 - (2*a_nonlinear*log(l_min/l_max) - (l_min - l_max))^2 - 2*log(l_min/l_max)*(l_min^2 - l_max^2)) / (4*log(l_min/l_max)*(l_max - l_min))
p_hookean = (k = ks, a = a_hookean)

F₁ = (δ,p) -> p.k .* (δ .- p.a) # Hookes law (Linear springs)
F₂ = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a .- (1 ./ δ)) # Nonlinear restoring force (Nonlinear springs)


hookean_area = Aₕ(l_min,l_max,p_hookean)
nonlinear_area = Aₙ(l_min,l_max,p_nonlinear)

hookean = F₁(l,p_hookean)
nonlinear = F₂(l,p_nonlinear)

function plotForceCompare(l, F_hookean, F_nonlinear, a_hookean, a_nonlinear)
    f = Figure(fontsize = 35,backgroundcolor=RGBf(1.0, 1.0, 1.0),
    resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    ymin = minimum([minimum(F_hookean), minimum(F_nonlinear)]) - 10
    ymax = maximum([maximum(F_hookean), maximum(F_nonlinear)]) + 10
    gaxmain = Axis(ga[1, 1], limits=(l[1], l[end], ymin, ymax), aspect=AxisAspect(1), 
                    xlabel=L"\text{Spring length [μm]}", ylabel=L"\text{Force Amplitude}")
    CairoMakie.lines!(gaxmain,l, F_hookean, linewidth=5, color=:blue)
    CairoMakie.lines!(gaxmain,l, F_nonlinear,linewidth=5, color=:green)
    CairoMakie.vlines!(gaxmain, [a_hookean, a_nonlinear], linewidth=5, linestyle=:dash, color=[:black, :grey])
    return f
end

f = plotForceCompare(l, hookean, nonlinear, a_hookean, a_nonlinear)

save("Force_Compare_2D_a_$a_nonlinear.png", f)


## relationship between a_nonlinear and a_hookean
a_h = (a_nonlinear, l_min, l_max) -> ((l_max - l_min)^2 .- (2 .*a_nonlinear*log(l_min/l_max) .- (l_min - l_max)).^2 .- 2 *log(l_min/l_max)*(l_min^2 - l_max^2)) ./ (4 *log(l_min/l_max)*(l_max - l_min))
nonlinear_resting_lengths_min = 5.0
nonlinear_resting_lengths_max = 20.0
nonlinear_resting_lengths = LinRange(nonlinear_resting_lengths_min, nonlinear_resting_lengths_max, 100)
linear_resting_lengths = a_h(nonlinear_resting_lengths,l_min,l_max)

function plotRestingLengthCompare(nonlinear_resting_lengths, linear_resting_lengths)
    f = Figure(fontsize = 35,backgroundcolor=RGBf(1.0, 1.0, 1.0),
    resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], aspect=AxisAspect(1), 
                    xlabel=L"\text{Nonlinear resting length [μm]}", ylabel=L"\text{Hookean resting length [μm]}")
    CairoMakie.lines!(gaxmain,nonlinear_resting_lengths, linear_resting_lengths, linewidth=5, color=:blue)
    return f
end

f_resting_lengths = plotRestingLengthCompare(nonlinear_resting_lengths, linear_resting_lengths)

## Stress Plots
ψ = (F, A) -> F./A

ψ_hookean = ψ(hookean,l)
ψ_nonlinear = ψ(nonlinear,l)


function plotStressCompare(l, ψ_hookean, ψ_nonlinear, a_hookean, a_nonlinear)
    f = Figure(fontsize = 35,backgroundcolor=RGBf(1.0, 1.0, 1.0),
    resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    ymin = minimum([minimum(ψ_hookean), minimum(ψ_hookean)]) - 10
    ymax = maximum([maximum(ψ_nonlinear), maximum(ψ_nonlinear)]) + 10
    gaxmain = Axis(ga[1, 1], limits=(l[1], l[end], ymin, ymax), aspect=AxisAspect(1), 
                    xlabel=L"\text{Spring length [μm]}", ylabel=L"\text{Stress [N/μm²]}")
    CairoMakie.lines!(gaxmain,l, ψ_hookean, linewidth=5, color=:blue)
    CairoMakie.lines!(gaxmain,l, ψ_nonlinear,linewidth=5, color=:green)
    CairoMakie.vlines!(gaxmain, [a_hookean, a_nonlinear], linewidth=5, linestyle=:dash, color=[:black, :grey])
    return f
end

f = plotStressCompare(l, ψ_hookean, ψ_nonlinear, a_hookean, a_nonlinear)

save("Stress_Compare_2D_a_$a0.png", f)





## Checking relationship between D and kₛ where D is a function of resting length a => D(a)

kₛ_func = (D,a,η) -> D.*η ./ (a.^2)

D = 0.0075
η = 1
a_min = 5
a_max =  15
a = LinRange(a_min, a_max, 100)

stiffness = kₛ_func(D,a,η)

function plotStiffnessVSRestingLength(a, stiffness)
    f = Figure(fontsize = 35,backgroundcolor=RGBf(1.0, 1.0, 1.0),
    resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    ymin = minimum(stiffness)
    ymax = maximum(stiffness)
    gaxmain = Axis(ga[1, 1], limits=(a[1], a[end], ymin, ymax), aspect=AxisAspect(1), 
                    xlabel=L"\text{Spring length [μm]}", ylabel=L"\text{Stiffness}",
                    title=L"D = 0.0075")
    
    CairoMakie.lines!(gaxmain, a, stiffness, linewidth=5, color=:blue)
    return f
end

f_stiffness = plotStiffnessVSRestingLength(a, stiffness)

