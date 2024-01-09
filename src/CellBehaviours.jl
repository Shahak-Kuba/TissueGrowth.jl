using LinearAlgebra

function calc_spring_densities(uᵢ)
    uᵢ₊₁ = circshift(uᵢ',1)
    return 1 ./ δ(uᵢ₊₁, uᵢ')
end

function calc_cell_densities(u,m)
    array = calc_spring_densities(u)
    return [sum(array[i:i+m-1]) for i in 1:m:length(array)-m+1]
end


P(event,ρ,α) = event ? α.*ρ : zeros(size(ρ))
A(event,ρ,β) = event ? β.*ρ : zeros(size(ρ))
E(event,ρ,γ) = event ? γ.*ρ : zeros(size(ρ))

function cell_probs(uᵢ,m,δt,prolif,death,embed,α,β,γ)
    ρ = calc_cell_densities(uᵢ,m)
    return (P(prolif,ρ,α).*δt, A(death, ρ,β).*δt, E(embed, ρ,γ).*δt)
end

function find_cell_index(arr::Vector{Float64}, threshold::Float64)
    cum_sum = cumsum(arr)
    index = findfirst(cum_sum .>= threshold)
    if index === nothing || index > length(arr)
        return index = length(arr)  # If the cumulative sum never reaches the threshold
    end
    return index
end



function affect!(integrator)
    (m,kₛ,η,kf,l₀,δt,growth_dir,prolif,death,embed,α,β,γ) = integrator.p
    u = integrator.u
    (p,a,e) = cell_probs(u, m, δt, prolif, death, embed, α, β, γ)
    (r1,r2,r3) = rand(3)
    if r1 < (sum(p) + sum(a) + sum(e)) # check if event has occurred
        if r2 < sum(p) / (sum(p) + sum(a) + sum(e)) # prolif occurs
            idx = find_cell_index(p, r3 * sum(p))
            println("inserted at: ",idx)
            insert!(u,idx,((circshift(u',1)'[:,idx] + u[:,idx])/2))
            # Perform operations based on prolif occurrence if needed
 

        elseif sum(p) / (sum(p) + sum(a) + sum(e)) < r2 < (sum(p) + sum(a)) / (sum(p) + sum(a) + sum(e)) # death occurs
            idx = find_cell_index(a, r3 * sum(a))
            println("deleted at: ",idx)
            deleteat!(u,idx)
            resize!(integrator,(2,size(integrator.u,2)))
            # Remove the cell from u based on death occurrence


        else # embed occurs
            idx = find_cell_index(e, r3 * sum(e))
            # Perform operations based on embed occurrence if needed
        end
        resize!(integrator,(2,size(integrator.u,2)))
    end
    nothing
end