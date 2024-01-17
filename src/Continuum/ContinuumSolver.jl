using LinearAlgebra
using Plots
using BenchmarkTools

include("SolverFncs.jl")

# simulation parameters
D = 1.0;
kf = 0.002
A = 0.0;
ρ₀ = 22.918311805232928;
growth_dir = "inward"
Tmax = 25.0
Xmax = 2π


function SolveContinuumLim1D(D,kf,A,ρ₀,growth_dir,Tmax,Xmax)
    # time and space discretisation
    T₀ = 0.0
    X₀ = 0.0

    N = 3000
    Δt = (Tmax - T₀)/N
    t = LinRange(T₀,Tmax,Int64((Tmax - T₀)/Δt))

    M = 200
    Δx = (Xmax - X₀)/M
    x = LinRange(X₀,Xmax,M)

    # initialise all derivatives and variables
    h = zeros(N+1,M)
    hₓ = zeros(N+1,M)
    hₓₓ = zeros(N+1,M)
    ρ = zeros(N+1,M)
    ρₓ = zeros(N+1,M)
    κ = zeros(N+1,M)
    λ = zeros(N+1,M)
    Φ = zeros(N+1,M)
    Ψ = zeros(N+1,M)

    # initialise matricies for cyclic tridiagonal
    TriDiagMat = zeros(M,M)
    d₀ = zeros(Float64,M)
    d₁ = zeros(Float64,M-1)
    d₋₁ = zeros(Float64,M-1)


    # set initial conditions for h and ρ
    h[1,:] = 2.0 .+ 0.5.*cos.(3*x)
    h[1,end] = h[1,1]
    #x = 3 .* cos.(x)
    ρ[1,:] = ones(1,M)*ρ₀

    if growth_dir == "inward"
        S = 1.0
    else
        S = -1.0
    end

    a⁺ = zeros(size(h[1,:]))
    a⁻ = zeros(size(h[1,:]))
    hᵢ = zeros(size(h[1,:]))
    hᵢ₊₁ = zeros(size(h[1,:]))
    hᵢ₋₁ = zeros(size(h[1,:]))
    ρᵢ = zeros(size(ρ[1,:]))
    ρᵢ₊₁ = zeros(size(ρ[1,:]))
    ρᵢ₋₁ = zeros(size(ρ[1,:]))
    
    @time for n in 1:N

        hᵢ .= h[n,:]
        hᵢ₊₁ .= circshift(hᵢ,1)
        hᵢ₋₁ .= circshift(hᵢ,-1)
        ρᵢ = ρ[n,:]
        ρᵢ₊₁ .= circshift(ρᵢ,1)
        ρᵢ₋₁ .= circshift(ρᵢ,-1)

        a⁺ .= aₘ₊(hᵢ₊₁,hᵢ₋₁)
        a⁻ .= aₘ₋(hᵢ₊₁,hᵢ₋₁)

        hₓ[n,:] .= (hᵢ.*(a⁺.-a⁻) .- hᵢ₋₁.*a⁺ .+ hᵢ₊₁.*a⁻) ./ Δx;                    # upwind to find hₓ,  refer eqn (26) in the notes
        hₓₓ[n,:] .= (hᵢ₊₁ .- 2 .*hᵢ .+ hᵢ₋₁) ./ (Δx^2);                               # cental to find hₓₓ, refer eqn (30) in the notes  
        ρₓ[n,:] .= (ρᵢ .* (a⁺.-a⁻) .- ρᵢ₋₁ .* a⁺ .+ ρᵢ₊₁ .* a⁻) ./ Δx;

        κ[n,:] = - hₓₓ[n,:] ./ ( (1 .+ hₓ[n,:].^2).^1.5 )

        Φ[n,:] = ρᵢ .- S.*Δt.*κ[n,:].*kf.*(ρᵢ.^2) .+ (S.*Δt.*kf.*ρᵢ.*ρₓ[n,:].*hₓ[n,:])./(sqrt.(1 .+(hₓ[n,:].^2))) .- D.*Δt.*ρₓ[n,:].*hₓ[n,:].*hₓₓ[n,:] ./ ((1 .+ (hₓ[n,:].^2)).^2) - A.*Δt.*ρᵢ
        λ[n,:] = (D.*Δt./(Δx.^2)) ./ (1 .+ hₓ[n,:].^2)

        h[n+1,:] = hᵢ .+ S.*Δt.*kf.*ρᵢ.*sqrt.(1 .+ (hₓ[n,:].^2))

        # cyclic tridiagonal matrix
        d₀ .= 1 .+ 2 .* λ[n,1:M]
        d₁ .= -λ[n,1:M-1]
        d₋₁ .= -λ[n,2:M]
    
        TriDiagMat .= diagm(-1 => d₋₁, 0 => d₀, 1 => d₁)
        TriDiagMat[M,1] = -λ[n,M]
        TriDiagMat[1,M] = -λ[n,1]

        # solving the system
        ρ[n+1,1:M] .= TriDiagMat\Φ[n,1:M]

    end

    return x,h,ρ

end



