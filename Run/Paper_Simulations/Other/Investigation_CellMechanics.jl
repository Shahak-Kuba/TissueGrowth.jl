using Plots

F₁ = (δ,p) -> p.k .* (δ .- p.a) # hookes law
F₂ = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a .- (1 ./ δ)) # nonlinear restoring force
F₂lim = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a)

x_min = 0.02
x_max = 0.034
x = LinRange(x_min, x_max, 100)

kₕ = 12.5
aₕ = 0.02
pₕ = (k = kₕ, a = aₕ)

hookean = F₁(x,pₕ)

kₙ = 12.5
aₙ = 0.02
pₙ = (k = kₙ, a = aₙ)

nonlinear = F₂(x,pₙ)
nonlinear_lim = F₂lim(x,pₙ)


#plotting
f = Plots.plot(x,[hookean nonlinear],linewidth=3,label=["Hookean" "Nonlinear"], xlabel="length [mm]", ylabel="Force [N/mm]")
Plots.plot!(x, nonlinear_lim,linewidth=3,linestyle=:dash, label="Nonlinear Limit")

save("Force_Compare.png", f)