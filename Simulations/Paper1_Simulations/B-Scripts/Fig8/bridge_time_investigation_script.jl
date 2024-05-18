using Makie
using CairoMakie

# Time to bridge based on side length
#Tbₛ = (sₛ, kf, q₀) -> sₛ./(4*kf*q₀)
#Tbₕ = (sₕ, kf, q₀) -> (√3 .* sₕ)./(4*kf*q₀)

#Ωmax = 500^2
#Ωmin = 50^2
#Ω₀ = LinRange(Ωmax, Ωmin, 100)
#Sₛ = sqrt.(Ω₀)
#Sₕ = sqrt.((2/(3√3)).*Ω₀)

#ratio = Sₛ ./ Sₕ

#q₀ = 1/20

#Tb_square = Tbₛ(Sₛ, kf, q₀)
#Tb_hex = Tbₕ(Sₕ, kf, q₀)

Tb_Kuba_2024 = (s, kf, q₀) -> s./(4*kf*q₀)
Tb_Buenzli_2020 = (s,Tb₀,μ) -> Tb₀.*s.^μ 


S = LinRange(100,700,500)

# Our analytical model

q₀ = 1/20;
kf = 87.842
Tb_2024 = Tb_Kuba_2024(S,kf,q₀)

# Buenzli et al. 2020
Tb₀_low = 0.050646 - 0.0056870
Tb₀_mid = 0.050646
Tb₀_high = 0.050646 + 0.0056870

μ_low = 0.999094 - 0.181650
μ_mid = 0.999094 
μ_high = 0.999094 + 0.181650

Tb_2020 = Tb_Buenzli_2020(S,Tb₀_mid, μ_high)

function plotAnalytic_vs_Regression(S, Tb_2020, Tb_2024)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 650, width=650, limits=(100,700,5,40),
                    xlabel="L [μm]", xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Tb [days]", ylabelsize = txtSize, yticklabelsize = tickSize)
    
    
    CairoMakie.lines!(gaxmain,S,Tb_2020,label="Buenzli et al. 2020", linewidth=4, linestyle=:dash, color=:black)
    CairoMakie.lines!(gaxmain,S,Tb_2024,label="Tb(S)", linewidth=4, linestyle=:solid, color=:blue)

    #Legend(f[1,1],[Analytic_Sol,Square_Sol,Hex_Sol], ["Analytic Circle", "Discrete Square","Discrete Hex"])
    axislegend(gaxmain, merge = true, unique = true, labelsize=tickSize)
    return f
end

f = plotAnalytic_vs_Regression(S, Tb_2020, Tb_2024)