KF = 8784.2;
Tb = 28.46
l = 500;
Ω₀ = l^2
P = l*4
q₀ = 1/20;
N = 100 #Int(P*q₀) # number of cells
kf = KF/N
l_min = 5
l_max = 20

# Time to bridge based on side length
Tbₛ = (sₛ, kf, q₀) -> sₛ./(4*kf*q₀)
Tbₕ = (sₕ, kf, q₀) -> (√3 .* sₕ)./(4*kf*q₀)

Ωmax = 500^2
Ωmin = 50^2
Ω₀ = LinRange(Ωmax, Ωmin, 100)
Sₛ = sqrt.(Ω₀)
Sₕ = sqrt.((2/(3√3)).*Ω₀)

ratio = Sₛ ./ Sₕ

q₀ = 1/20

Tb_square = Tbₛ(Sₛ, kf, q₀)
Tb_hex = Tbₕ(Sₕ, kf, q₀)