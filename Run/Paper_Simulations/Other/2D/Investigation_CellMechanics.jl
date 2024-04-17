using Plots

F₁ = (δ,p) -> p.k .* (δ .- p.a) # hookes law
F₂ = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a .- (1 ./ δ)) # nonlinear restoring force
F₂lim = (δ,p) -> p.k .* p.a^2 .* (ones(size(δ))./p.a)

x_min = 4
x_max = 21
x = LinRange(x_min, x_max, 100)

kₕ = 1.25e6
aₕ = 10
pₕ = (k = kₕ, a = aₕ)

hookean = F₁(x,pₕ)

kₙ = 1.25e6
aₙ = 10
pₙ = (k = kₙ, a = aₙ)

nonlinear = F₂(x,pₙ)
nonlinear_lim = F₂lim(x,pₙ)


#plotting
f = Plots.plot(x,[hookean nonlinear],linewidth=3,label=["Hookean" "Nonlinear"], xlabel="cell length [μm]", ylabel="Force Amplitude")
Plots.vline!([10.0], linewidth=3, linestyle=:dash,label="Resting length")
Plots.vline!([5], linewidth=3, linestyle=:dash,label="Minimum length")
Plots.vline!([20.0], linewidth=3, linestyle=:dash,label="Maximum length")
#Plots.plot!(x, nonlinear_lim,linewidth=3,linestyle=:dash, label="Nonlinear Limit")

save("Force_Compare_2D_Infil.png", f)