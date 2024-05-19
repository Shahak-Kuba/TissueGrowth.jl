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

Tb_Square = (s, kf, q₀) -> s./(4*kf*q₀)
Tb_Hex = (s, kf, q₀) -> (√3 .*s)./(4*kf*q₀)
Tb_Buenzli_2020 = (s,Tb₀,μ) -> Tb₀.*s.^μ 


S = LinRange(100,700,500)

# Our analytical model

q₀ = 1/20;
kf = 87.842
Tb_Square_2024 = Tb_Square(S,kf,q₀)
Tb_Hex_2024 = Tb_Hex(S,kf,q₀)

# Buenzli et al. 2020
Tb₀_low = 0.050646 - 0.0056870
Tb₀_mid = 0.050646
Tb₀_high = 0.050646 + 0.0056870

μ_low = 0.999094 - 0.181650
μ_mid = 0.999094 
μ_high = 0.999094 + 0.181650

Tb_2020 = Tb_Buenzli_2020(S,Tb₀_high, μ_mid)  # Nice Fit
#Tb_2020 = Tb_Buenzli_2020(S,Tb₀_mid, μ_mid)

function plotAnalytic_vs_Regression(S, Tb_2020, Tb_Square_2024, Tb_Hex_2024)
    txtSize = 35;
    tickSize = 25;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()
    gaxmain = Axis(ga[1, 1], height = 650, width=650, limits=(100,700,0,75),
                    xlabel="L [μm]", xlabelsize = txtSize, xticklabelsize = tickSize,
                    ylabel="Tb [days]", ylabelsize = txtSize, yticklabelsize = tickSize)
    
    
    CairoMakie.lines!(gaxmain,S,Tb_Square_2024,label=L"\text{Square:}\;T_{b}(L)", linewidth=4, linestyle=:solid, color=:blue)
    CairoMakie.lines!(gaxmain,S,Tb_Hex_2024,label=L"\text{Hex:}\;T_{b}(L)", linewidth=4, linestyle=:solid, color=:red)
    CairoMakie.lines!(gaxmain,S,Tb_2020,label=L"\text{Buenzli et al. 2020}", linewidth=4, linestyle=:dash, color=:black)


    #Legend(f[1,1],[Analytic_Sol,Square_Sol,Hex_Sol], ["Analytic Circle", "Discrete Square","Discrete Hex"])
    axislegend(gaxmain, merge = true, unique = true, labelsize=tickSize, position=:lt)
    return f
end

f = plotAnalytic_vs_Regression(S, Tb_2020, Tb_Square_2024, Tb_Hex_2024)