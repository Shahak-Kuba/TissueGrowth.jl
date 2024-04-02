# Analytic solution for number of cells
N_analytic(N₀,A,kf,Ot,t) = N₀.*exp.(-(A+kf*Ot).*t)

# Analytic solution for void area
Ω_analytic(Ω₀,N₀,kf,t) = Ω₀ .- kf.*N₀.*t # constant cell population
Ω_analytic(Ω₀,N₀,A,kf,Ot,t) = Ω₀ .- (N₀./(A.+kf.*Ot)).*(1 .- exp.(-(A+kf*Ot).*t)) # changing cell population

# function to calculate kf value from experiment

Rr(KF, Ω₀, t)= sqrt.((Ω₀.-KF.*t)/π)

V(KF, Ω₀, t) = (-KF ./ (2 .*(π^(0.5))) .* ((1)./((Ω₀ .- KF.*t).^(0.5))))

V(KF, Ω₀, Tb, t) = (KF ./ 2π) .* ((Ω₀ - KF.*(Tb .- t))./(π)).^(0.5)

Ncells(KF, Ω₀, ρ₀, t) = 2 .* π .* ρ₀ .* ((Ω₀ .- KF.*t) ./ π).^(0.5)