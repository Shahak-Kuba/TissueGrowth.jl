using TissueGrowth
using Makie
using CairoMakie

## Resting length investigation (with Phase shift)
#Φ = (l₀, l₁, p) -> (p.k * ( (2*p.a*(1/l₀ - 1/l₁)) - ((p.a^2 *((1/l₀^2 - 1/l₁^2)))/2) + log((1/l₁)/(1/l₀)) ) ) / (1/l₁ - 1/l₀)
#Aₙ = (l₀, l₁, ξ, p) -> p.k * p.a^2 * ((1/l₁ - 1/l₀)/p.a + (1/l₀^2 - 1/l₁^2)/2) + ξ*(1/l₁ - 1/l₀)
#F₂ = (δ,p,ξ) -> p.k .* p.a^2 .* (ones(size(δ))./p.a .- (1 ./ δ)) .+ ξ # nonlinear restoring force w/ verticle shift
#Φ_value = Φ(l_min,l_max,p)
#a1 = 1 / (1/a0 + Φ_value/(ks*a0^2))
#a1_2 = ks*a0^2 / (ks*a0 + Φ_value)

## Resting length investigation (with a₁ and a₂ + no phase shift)
Aₕ = (l₀, l₁, p) -> p.k .* ((l₁.^2 - l₀.^2)./2 + p.a .* (l₀ - l₁))
Aₙ = (l₀, l₁, p) -> p.k .* p.a^2 .* ((l₁ .- l₀)./p.a .+ log.(l₀) .- log.(l₁))

l_min = 2.5
l_max = 20.0
l = LinRange(l_min, l_max, 100)

ks = 20.0
a0 = 10.0
p = (k = ks, a = a0)

a1 = (l_min - l_max - √((l_max - l_min)^2 - 4*log(l_min/l_max)*(a0*(l_max - l_min) + (l_min^2 - l_max^2)/2))) / (2*log(l_min/l_max))
pₙ = (k = ks, a = a1)

F₁ = (δ,p) -> p.k .* (δ .- p.a) # hookes law
F₂ = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a .- (1 ./ δ))


hookean_area = Aₕ(l_min,l_max,p)
nonlinear_area = Aₙ(l_min,l_max,pₙ)

hookean = F₁(l,p)
nonlinear = F₂(l,pₙ)

function plotForceCompare(l, F_hookean, F_nonlinear, a0)
    f = Figure(fontsize = 35,backgroundcolor=RGBf(1.0, 1.0, 1.0),
    resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    ymin = minimum([minimum(F_hookean), minimum(F_nonlinear)]) - 10
    ymax = maximum([maximum(F_hookean), maximum(F_nonlinear)]) + 10
    gaxmain = Axis(ga[1, 1], limits=(l[1], l[end], ymin, ymax), aspect=AxisAspect(1), 
                    xlabel=L"\text{Spring length [μm]}", ylabel=L"\text{Force Amplitude}")
    CairoMakie.lines!(gaxmain,l, F_hookean, linewidth=5, color=:blue)
    CairoMakie.lines!(gaxmain,l, F_nonlinear,linewidth=5, color=:green)
    CairoMakie.vlines!(gaxmain, [a0, a1], linewidth=5, linestyle=:dash, color=[:black, :grey])
    return f
end

f = plotForceCompare(l, hookean, nonlinear, a0)

save("Force_Compare_2D_a_$a0.png", f)



## Stress Plots
ψ = (F, A) -> F./A

ψ_hookean = ψ(hookean,l)
ψ_nonlinear = ψ(nonlinear,l)


function plotStressCompare(l, ψ_hookean, ψ_nonlinear, a0)
    f = Figure(fontsize = 35,backgroundcolor=RGBf(1.0, 1.0, 1.0),
    resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    ymin = minimum([minimum(ψ_hookean), minimum(ψ_hookean)]) - 10
    ymax = maximum([maximum(ψ_nonlinear), maximum(ψ_nonlinear)]) + 10
    gaxmain = Axis(ga[1, 1], limits=(l[1], l[end], ymin, ymax), aspect=AxisAspect(1), 
                    xlabel=L"\text{Spring length [μm]}", ylabel=L"\text{Stress [N/μm²]}")
    CairoMakie.lines!(gaxmain,l, ψ_hookean, linewidth=5, color=:blue)
    CairoMakie.lines!(gaxmain,l, ψ_nonlinear,linewidth=5, color=:green)
    CairoMakie.vlines!(gaxmain, [a0, a1], linewidth=5, linestyle=:dash, color=[:black, :grey])
    return f
end

f = plotStressCompare(l, ψ_hookean, ψ_nonlinear, a0)

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

