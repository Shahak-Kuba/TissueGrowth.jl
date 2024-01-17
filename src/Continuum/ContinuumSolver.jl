using LinearAlgebra

include("SolverFncs.jl")

# simulation parameters
D = 0.02;
A = 0.0;
growth_dir = "inward"

# time and space discretisation
T₀ = 0.0
Tmax = 250.0
X₀ = 0.0
Xmax = 2π

N = 3000
Δt = (Tmax - T₀)/N
t = LinRange(T₀,Tmax,Int64((Tmax - T₀)/Δt))

M = 200
Δx = (Xmax - X₀)/M
x = LinRange(X₀,Xmax,Int64(round((Xmax - X₀)/Δx)))

# initialise all derivatives and variables
h = zeros(N+1,M+1)
hₓ = zeros(N+1,M+1)
hₓₓ = zeros(N+1,M+1)
ρ = zeros(N+1,M+1)
ρₓ = zeros(N+1,M+1)
κ = zeros(N+1,M+1)
λ = zeros(N+1,M+1)
Φ = zeros(N+1,M+1)
Ψ = zeros(N+1,M+1)

# initialise matricies for cyclic tridiagonal
TriDiagMat = zeros(M,M)
d₀ = zeros(1,M)
d₁ = zeros(1,M-1)
d₋₁ = zeros(1,M-1)


# set initial conditions for h and ρ
h[1,1:end-1] = 2.0 .+ 0.5.*cos.(3*x)
h[1,end] = h[1,1]
ρ[1,:] = ones(1,M+1)*0.016

if growth_dir == "inward"
    S = 1.0
else
    S = -1.0
end
H = 1.0
a⁺ = 0.0
a⁻ = 0.0

for n in 1:1
    for m in 1:M
        if m == 1
            a⁺ = aₘ₊(h[n,2],h[n,M])
            a⁻ = aₘ₋(h[n,2],h[n,M])
            hₓ[n,1] = (h[n,1]*(a⁺-a⁻) - h[n,M]*a⁺ + h[n,2]*a⁻)/Δx;                       # upwind to find hₓ,  refer eqn (26) in the notes
            hₓₓ[n,1] = (h[n,2]-2*h[n,1]+h[n,M])/(Δx^2);                                  # cental to find hₓₓ, refer eqn (30) in the notes  
            ρₓ[n,1] = (ρ[n,1]*(a⁺-a⁻) - ρ[n,M]*a⁺ + ρ[n,2]*a⁻)/Δx;                       # upwind to find ρₓ,  refer eqn (29) in the notes  
        else 
            a⁺ = aₘ₊(h[n,m+1],h[n,m-1])
            a⁻ = aₘ₋(h[n,m+1],h[n,m-1])
            hₓ[n,m] = (h[n,m]*(a⁺-a⁻) - h[n,m-1]*a⁺ + h[n,m+1]*a⁻)/Δx;                    # upwind to find hₓ,  refer eqn (26) in the notes
            hₓₓ[n,m] = (h[n,m+1]-2*h[n,m]+h[n,m-1])/(Δx^2);                               # cental to find hₓₓ, refer eqn (30) in the notes  
            ρₓ[n,m] = (ρ[n,m]*(a⁺-a⁻) - ρ[n,m-1]*a⁺ + ρ[n,m+1]*a⁻)/Δx;                    # upwind to find ρₓ,  refer eqn (29) in the notes  
        end

        κ[n,m] = - hₓₓ[n,m] / ( (1 + hₓ[n,m]^2)^1.5 )

        # can apply heaviside function in case want only κ ≥ 0 to move in the longtitudinal direction (For H variable)

        Φ[n,m] = ρ[n,m] - S*Δt*κ[n,m]*H*(ρ[n,m]^2) + (S*Δt*H*ρ[n,m]*ρₓ[n,m]*hₓ[n,m])/(sqrt(1+(hₓ[n,m]^2))) - D*Δt*ρₓ[n,m]*hₓ[n,m]*hₓₓ[n,m] / ((1 + (hₓ[n,m]^2))^2) - A*Δt*ρ[n,m]
        λ[n,m] = (D*Δt/(Δx^2)) / (1 + hₓ[n,m]^2)

        h[n+1,m] = h[n,m] + S*Δt*H*ρ[n,m]*sqrt(1 + (hₓ[n,m]^2))
    end

    # cyclic tridiagonal matrix
    d₀ = 1 .+ 2 .* λ[n,1:M]
    d₁ = -λ[n,1:M-1]
    d₋₁ = -λ[n,2:M]
   
    TriDiagMat = diagm(-1 => d₋₁, 0 => d₀, 1 => d₁)
    TriDiagMat[M,1] = -λ[n,M]
    TriDiagMat[1,M] = -λ[n,1]

    # solving the system
    ρ[n+1,1:M] = TriDiagMat\Φ[n,1:M]

    # periodic boundary conditions
    h[n+1,M+1] = h[n+1,1]
    h[n+1,M+1] = ρ[n+1,1]
end

plot(x,h[1,1:end-1])
plot!(x,h[10,1:end-1])
plot!(x,h[20,1:end-1])
plot!(x,h[30,1:end-1])