using Plots

F₁ = (δ,p) -> p.k .* (δ .- p.a) # hookes law
F₂ = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a .- (1 ./ δ)) # nonlinear restoring force
F₂lim = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a)

kₙ = 0.005
aₙ = 1.0
pₙ = (k = kₙ, a = aₙ)

kₕ = 0.1
aₕ = 1.0
pₕ = (k = kₕ, a = aₕ)

x_min = 0.01
x_max = 0.05
x = LinRange(x_min, x_max, 100)

hookean = F₁(x,pₕ)
nonlinear = F₂(x,pₙ)
nonlinear_lim = F₂lim(x,p)

#plotting
Plots.plot(x,[hookean nonlinear],linewidth=3,label=["Hookean" "Nonlinear"])
Plots.plot!(x, nonlinear_lim,linewidth=3,linestyle=:dash, label="Nonlinear Limit")