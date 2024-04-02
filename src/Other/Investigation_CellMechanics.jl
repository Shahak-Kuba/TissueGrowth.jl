using Plots

F₁ = (δ,p) -> p.k .* (δ .- p.a) # hookes law
F₂ = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a .- (1 ./ δ)) # nonlinear restoring force
F₂lim = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a)

kₛ = 10.0
aₛ = 1.0
p = (k = kₛ, a = aₛ)

x_min = 0.25
x_max = 5.0
x = LinRange(x_min, x_max, 100)

hookean = F₁(x,p)
nonlinear = F₂(x,p)
nonlinear_lim = F₂lim(x,p)

#plotting
Plots.plot(x,[hookean nonlinear],linewidth=3,label=["Hookean" "Nonlinear"])
Plots.plot!(x, nonlinear_lim,linewidth=3,linestyle=:dash, label="Nonlinear Limit")