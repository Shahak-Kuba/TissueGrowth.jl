using Plots

include("ContinuumSolver.jl")

# simulation parameters
D = 0.0001;
kf = 0.002
A = 0.0;
ρ₀ = 22.918311805232928;
growth_dir = "inward"
Tmax = 25.0
Xmax = 2π

x,h,ρ = SolveContinuumLim1D(D,kf,A,ρ₀,growth_dir,Tmax,Xmax);

θ,R,ρ = SolveContinuumLim_Polar(D,kf,A,ρ₀,growth_dir,Tmax)

plot(x,h[1,:])
plot!(x,h[500,:])
plot!(x,h[1000,:])
plot!(x,h[1500,:])
plot!(x,h[2000,:])
plot!(x,h[2500,:])

plot(R[1,:].*cos.(θ),R[1,:].*sin.(θ))
plot!(R[500,:].*cos.(θ),R[500,:].*sin.(θ))
plot!(R[2500,:].*cos.(θ),R[2500,:].*sin.(θ))

plot(x,ρ[1,1:end])
plot!(x,ρ[500,1:end])
plot!(x,ρ[1000,1:end])
plot!(x,ρ[1500,1:end])
plot!(x,ρ[2000,1:end])
plot!(x,ρ[2500,1:end])