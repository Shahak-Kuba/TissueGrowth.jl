using Plots

## Resting length investigation
Φ = (l₀, l₁, p) -> (p.k * ( (2*p.a*(1/l₀ - 1/l₁)) - ((p.a^2 *((1/l₀^2 - 1/l₁^2)))/2) + log((1/l₁)/(1/l₀)) ) ) / (1/l₁ - 1/l₀)

Aₕ = (l₀, l₁, p) -> p.k * (log((1/l₁)/(1/l₀)) + p.a*((1/l₀) - (1/l₁)))
Aₙ = (l₀, l₁, ξ, p) -> p.k * p.a^2 * ((1/l₁ - 1/l₀)/p.a + (1/l₀^2 - 1/l₁^2)/2) + ξ*(1/l₁ - 1/l₀)

l_min = 7.5
l_max = 32.5
l = LinRange(l_min, l_max, 100)

ks = 12.5
a0 = 20
p = (k = ks, a = a0)

F₁ = (δ,p) -> p.k .* (δ .- p.a) # hookes law

Φ_value = Φ(l_min,l_max,p)
F₂ = (δ,p,ξ) -> p.k .* p.a^2 .* (ones(size(δ))./p.a .- (1 ./ δ)) .+ ξ # nonlinear restoring force w/ verticle shift

hookean_area = Aₕ(l_min,l_max,p)
nonlinear_area = Aₙ(l_min,l_max,Φ_value,p)

hookean = F₁(l,p)
nonlinear = F₂(l,p,Φ_value)


#plotting
f = Plots.plot(l,[hookean nonlinear],linewidth=3,label=["Hookean" "Nonlinear"], xlabel="cell length [μm]", ylabel="Force Amplitude")
Plots.vline!([20], linewidth=3, linestyle=:dash,label="Resting length")
Plots.vline!([10], linewidth=3, linestyle=:dash,label="Minimum length")
Plots.vline!([30], linewidth=3, linestyle=:dash,label="Maximum length")
#Plots.plot!(x, nonlinear_lim,linewidth=3,linestyle=:dash, label="Nonlinear Limit")

save("Force_Compare_1D_Trench.png", f)