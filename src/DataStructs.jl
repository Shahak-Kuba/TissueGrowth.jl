
using ElasticArrays


"""
SimResults_t(btype,t,u,∑F, ρ, V, Ω)

"""

struct SimResults_t
    btype::String
    t::Vector{Float64}
    u::Vector{ElasticMatrix{Float64,Vector{Float64}}}
    ∑F::Vector{ElasticVector{Float64,Vector{Float64}}}
    Density::Vector{ElasticMatrix{Float64,Vector{Float64}}}
    Vₙ::Vector{Vector{Float64}}
    Ω::Vector{Float64}
    ψ::Vector{ElasticMatrix{Float64,Vector{Float64}}}
    Κ::Vector{ElasticMatrix{Float64,Vector{Float64}}}
end