using Plots

l_min = 5
l_max = 20
kₛ = 10
Kₛ = 100
l₀ =20
ll = 0:0.1:l_max

#L = (l_min,l_max,kₛ,Kₛ,l₀) -> ((l_max .- l_min)./((kₛ./Kₛ).*((l_max^2 .- l_min^2)/2 .+ l₀.*(l_min .- l_max)) .- log.(l_min./l_max)))
L = (l_min,l_max,kₛ,Kₛ,l₀) -> (l_max .- l_min) ./ ((kₛ./Kₛ).*((l_max.^2 - l_min.^2)./2 .- l₀.*l_max .+ l₀.*l_min) .+ log.(l_max) .- log.(l_min))

L(l_min,l_max,kₛ,Kₛ,l₀)

f = Plots.plot(ll,L(l_min,5,kₛ,Kₛ,ll),linewidth=3,label="l_max=5")
for ii in 10:5:l_max
   Plots.plot!(f,ll,L(l_min,ii,kₛ,Kₛ,ll),linewidth=3,label="l_max=$ii")
end
f

#Kₛ = (l_min,l_max,kₛ,a) -> ( kₛ.*( (l_max.^2 - l_min.^2)./2) .- a.*(l_max-l_min) ) ./ ( (l_max - l_min)./a .+ log.(l_min./l_max))

#Kₛ = (l_min,l_max,kₛ,a) -> ( kₛ.*((l_max.^2 - l_min.^2)./2 .- a.*l_max .+ a.*l_min) ) ./ ( ((l_max .- l_min)./a .- log.(l_max) .+ log.(l_min)))
#Kₛ(l_min,l_max,kₛ,a)