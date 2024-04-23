using Plots

## Resting length investigation
Φ = (l₀, l₁, p) -> (p.k * ( (2*p.a*(1/l₀ - 1/l₁)) - ((p.a^2 *((1/l₀^2 - 1/l₁^2)))/2) + log((1/l₁)/(1/l₀)) ) ) / (1/l₁ - 1/l₀)

Aₕ = (l₀, l₁, p) -> p.k * (log((1/l₁)/(1/l₀)) + p.a*((1/l₀) - (1/l₁)))
Aₙ = (l₀, l₁, ξ, p) -> p.k * p.a^2 * ((1/l₁ - 1/l₀)/p.a + (1/l₀^2 - 1/l₁^2)/2) + ξ*(1/l₁ - 1/l₀)

l_min = 5.0
l_max = 20.0
l = LinRange(l_min, l_max, 100)

ks = 1.25
a0 = 10.0
p = (k = ks, a = a0)

F₁ = (δ,p) -> p.k .* (δ .- p.a) # hookes law

Φ_value = Φ(l_min,l_max,p)
F₂ = (δ,p,ξ) -> p.k .* p.a^2 .* (ones(size(δ))./p.a .- (1 ./ δ)) .+ ξ # nonlinear restoring force w/ verticle shift

hookean_area = Aₕ(l_min,l_max,p)
nonlinear_area = Aₙ(l_min,l_max,Φ_value,p)

hookean = F₁(l,p)
nonlinear = F₂(l,p,Φ_value)


#plotting
f = Plots.plot(1 ./l,[hookean nonlinear],linewidth=3,label=["Hookean" "Nonlinear"], xlabel="cell density [1/μm]", ylabel="Force Amplitude")
Plots.vline!([1/a0], linewidth=3, linestyle=:dash,label="Resting length")
Plots.vline!([1/5], linewidth=3, linestyle=:dash,label="Minimum length")
Plots.vline!([1/20], linewidth=3, linestyle=:dash,label="Maximum length")
#Plots.plot!(x, nonlinear_lim,linewidth=3,linestyle=:dash, label="Nonlinear Limit")

f = plotForceCompare(l, hookean, nonlinear, a0)

save("Force_Compare_2D_Square_a_$a0.png", f)


function plotForceCompare(l, F_hookean, F_nonlinear, a0)
    f = Figure(fontsize = 32,backgroundcolor=RGBf(1.0, 1.0, 1.0),
    resolution=(1000, 800))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], limits=(4, l[end] + 1, -10, 13), xlabel="Spring length [μm]", ylabel="Force Amplitude")
    CairoMakie.lines!(gaxmain,l, F_hookean, linewidth=5, color=:red)
    CairoMakie.lines!(gaxmain,l, F_nonlinear,linewidth=5)
    CairoMakie.vlines!(gaxmain, [l[1], a0, l[end]], linewidth=5, linestyle=:dash, color=[:grey, :black, :grey])
    return f
end