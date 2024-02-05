

function Set_Random_Seed(seednum=123)
    Random.seed!(seednum)
end

"""
    calc_spring_densities(uᵢ)

Calculate the spring densities for a given state vector `uᵢ`.

This function computes the inverse of the distance between each pair of adjacent elements in `uᵢ`. The calculation involves shifting `uᵢ` and then applying the δ function to each pair.

# Arguments
- `uᵢ`: A state vector representing the positions of particles or cells.

# Returns
A vector of spring densities, where each element is the inverse of the distance between adjacent elements in `uᵢ`.
"""
function calc_spring_densities(uᵢ)
    uᵢ₊₁ = circshift(uᵢ',1)
    return 1 ./ δ(uᵢ₊₁, uᵢ')
end

"""
    calc_cell_densities(u, m)

Calculate the cell densities over a window of size `m` for a given state vector `u`.

This function first computes the spring densities using `calc_spring_densities` and then calculates the sum of these densities over a sliding window of size `m`. 

# Arguments
- `u`: A state vector representing the positions of particles or cells.
- `m`: The size of the window over which to calculate the densities.

# Returns
A vector of cell densities, each computed over a window of size `m`.
"""
function calc_cell_densities(u,m)
    array = calc_spring_densities(u)
    return [sum(array[i:i+m-1]) for i in 1:m:length(array)-m+1]
end

"""
    P(event, ρ, α)

Calculate the probability of proliferation occurring, given the density `ρ` and the parameter `α`.

If the proliferation event is considered to occur (`event` is `true`), the probability is calculated as the product of `ρ` and `α`. Otherwise, a vector of zeros is returned.

# Arguments
- `event`: A boolean indicating whether the proliferation event is considered to occur.
- `ρ`: A vector representing densities.
- `α`: Parameter for scaling the proliferation probability.

# Returns
A vector representing the calculated probabilities for proliferation.
"""
P(event,ρ,α) = event ? α.*ρ : zeros(size(ρ))

"""
    A(event, ρ, β)

Calculate the probability of apoptosis (cell death) occurring, given the density `ρ` and the parameter `β`.

If the apoptosis event is considered to occur (`event` is `true`), the probability is calculated as the product of `ρ` and `β`. Otherwise, a vector of zeros is returned.

# Arguments
- `event`: A boolean indicating whether the apoptosis event is considered to occur.
- `ρ`: A vector representing densities.
- `β`: Parameter for scaling the apoptosis probability.

# Returns
A vector representing the calculated probabilities for apoptosis.
"""
A(event,ρ,β) = event ? β.*ρ : zeros(size(ρ))

"""
    E(event, ρ, γ)

Calculate the probability of embedding occurring, given the density `ρ` and the parameter `γ`.

If the embedding event is considered to occur (`event` is `true`), the probability is calculated as the product of `ρ` and `γ`. Otherwise, a vector of zeros is returned.

# Arguments
- `event`: A boolean indicating whether the embedding event is considered to occur.
- `ρ`: A vector representing densities.
- `γ`: Parameter for scaling the embedding probability.

# Returns
A vector representing the calculated probabilities for embedding.
"""
E(event,ρ,kf,γ) = event ? γ.*kf.*ones(size(ρ)) : zeros(size(ρ))

"""
    cell_probs(uᵢ, m, δt, prolif, death, embed, α, β, γ)

Calculate the probabilities for cell-related events (proliferation, death, embedding) based on the current state vector `uᵢ`.

This function uses `calc_cell_densities` to calculate cell densities and then computes the probabilities for proliferation, death, and embedding events.

# Arguments
- `uᵢ`: A state vector representing the positions of particles or cells.
- `m`: The size of the window over which to calculate the densities.
- `δt`: Time step size.
- `prolif`, `death`, `embed`: Boolean flags indicating whether the respective event is considered to occur.
- `α, β, γ`: Parameters for scaling the probabilities.

# Returns
A tuple containing three vectors representing the probabilities for proliferation, death, and embedding.
"""
function cell_probs(uᵢ,m,δt,prolif,death,embed,α,β,γ,kf)
    ρ = calc_cell_densities(uᵢ,m)
    return (P(prolif,ρ,α).*δt, A(death,ρ,β).*δt, E(embed,ρ,kf,γ).*δt)
end

"""
    find_cell_index(arr::Vector{Float64}, threshold::Float64)

Find the index in `arr` where the cumulative sum first exceeds or equals `threshold`.

This function is typically used in stochastic processes to determine an outcome based on a probability distribution.

# Arguments
- `arr`: A vector of probabilities.
- `threshold`: A threshold value used to find the corresponding index in `arr`.

# Returns
The index in `arr` where the cumulative sum first exceeds or equals `threshold`. If the threshold is not met, the length of `arr` is returned.

"""
function find_cell_index(arr::Vector{Float64}, threshold::Float64)
    cum_sum = cumsum(arr)
    index = findfirst(cum_sum .>= threshold)
    if index === nothing || index > length(arr)
        return index = length(arr)  # If the cumulative sum never reaches the threshold
    end
    return index
end

"""
    affect!(integrator)

Update the state of `integrator` based on probabilistic cellular events.

This function modifies `integrator` in place. It uses the parameters and state from `integrator` to compute probabilities for different cellular events: proliferation (`prolif`), death (`death`), and embedding (`embed`). Based on these probabilities and random draws, it updates the state of `integrator`.

# Arguments
- `integrator`: The integrator object containing the current state and parameters. The parameters are expected to be a tuple containing:

# Returns
`nothing`. The function modifies `integrator` in place.
"""
function affect!(integrator)
    (m,kₛ,η,kf,l₀,δt,growth_dir,prolif,death,embed,α,β,γ) = integrator.p
    u = integrator.u
    (p,a,e) = cell_probs(u, m, δt, prolif, death, embed, α, β, γ, kf)
    (r1,r2,r3) = rand(3)
    if r1 < (sum(p) + sum(a) + sum(e)) # check if event has occurred
        if r2 < sum(p) / (sum(p) + sum(a) + sum(e)) # prolif occurs
            idx = find_cell_index(p, r3 * sum(p))
            #println("inserted at: ",idx)
            insert!(u,idx,((circshift(u',1)'[:,idx] + u[:,idx])/2))
            # Perform operations based on prolif occurrence if needed
 

        elseif sum(p) / (sum(p) + sum(a) + sum(e)) < r2 < (sum(p) + sum(a)) / (sum(p) + sum(a) + sum(e)) # death occurs
            idx = find_cell_index(a, r3 * sum(a))
            #println("deleted at: ",idx)
            deleteat!(u,idx)
            resize!(integrator,(2,size(integrator.u,2)))
            # Remove the cell from u based on death occurrence


        else # embed occurs
            idx = find_cell_index(e, r3 * sum(e))
            #println("embedded at: ",idx)
            store_embedded_cell(u, idx, m)

            for i = idx:(idx + m)-1
                deleteat!(u,idx)
            end
            resize!(integrator,(2,size(integrator.u,2)))
            # Perform operations based on embed occurrence if needed
        end
        resize!(integrator,(2,size(integrator.u,2)))
    end
    nothing
end

"""
    store_embed_cell_pos(pos)

Stores the position of embedded cells into embedded_cells array.

This function inserts the position vector `pos` of a boundary of an embedded cell into the embedded_cells array.

# Arguments
- `pos`: position vector of the cell that is being embedded into the tissue should contain `m+1` values given `m` springs in the simulation

# Returns
`nothing`. The function modifies `embedded_cells` in place.
"""
function store_embed_cell_pos(pos)
    global embedded_cells
    insert!(embedded_cells,size(embedded_cells,2),pos)
    return nothing
end

function store_embedded_cell(u, idx, m)
    for i = idx:idx+m
        store_embed_cell_pos(u[:,idx].data)
    end
    return nothing
end

function store_embedded_cell_count(u, t, integrator)
    (m,kₛ,η,kf,l₀,δt,growth_dir,prolif,death,embed,α,β,γ) = integrator.p
    global embedded_cells
    #global embedded_cell_count
    cell_count = size(hcat(embedded_cells...),2)/(m+1)
    #push!(embedded_cell_count,cell_count)
    return cell_count
end