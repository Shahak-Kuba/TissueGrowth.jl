# Analytic solution for number of cells
N_analytic(N₀,A,kf,Ot,t) = N₀.*exp.(-(A+kf*Ot).*t)

# Analytic solution for void area
Ω_analytic(Ω₀,N₀,kf,t) = Ω₀ .- kf.*N₀.*t # constant cell population
Ω_analytic(Ω₀,N₀,A,kf,Ot,t) = Ω₀ .- (N₀./(A.+kf.*Ot)).*(1 .- exp.(-(A+kf*Ot).*t)) # changing cell population