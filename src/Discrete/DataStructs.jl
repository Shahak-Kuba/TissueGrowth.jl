
"""
    SimResults_t

A data structure representing the results of a simulation which can then be used for visualisation with the plotting code `PlottingFncs1D.jl` and `PlottingFncs2D.jl`.

This structure contains various fields that store the results and relevant data from a simulation process. It includes the type of simulation, time steps, state vectors, forces, densities, and other relevant quantities.

# Fields
- `btype`: A `String` indicating the type of simulation or boundary condition.
- `t`: A `Vector{Float64}` representing the time steps at which the simulation results are recorded.
- `u`: A `Vector` of `ElasticMatrix{Float64,Vector{Float64}}` representing the state of the system at each time step.
- `∑F`: A `Vector` of `ElasticVector{Float64,Vector{Float64}}`, representing the sum of forces at each time step.
- `Density`: A `Vector` of `ElasticMatrix{Float64,Vector{Float64}}` representing the density of the system at each time step.
- `Vₙ`: A `Vector{Vector{Float64}}` representing the normal velocity at each time step.
- `Ω`: A `Vector{Float64}` representing the void area within the pore at each time step.
- `ψ`: A `Vector` of `ElasticMatrix{Float64,Vector{Float64}}`, representing the stress of the springs/ cells at each time step.
- `Κ`: A `Vector` of `ElasticMatrix{Float64,Vector{Float64}}`, representing the approximate curvature of the moving boundary at each time step.
- `CellCount`: A `Vector{Int64}` representing the amount of active cells at a give time `t`.

"""
struct SimResults_t
    btype::String
    t::Vector{Float64}
    u::Vector{ElasticMatrix{Float64,Vector{Float64}}}
    ∑F::Vector{ElasticVector{Float64,Vector{Float64}}}
    Density::Vector{ElasticMatrix{Float64,Vector{Float64}}}
    Vₙ::Vector{Vector{Float64}}
    Ω::Vector{Float64}
    ψ::Vector{ElasticVector{Float64,Vector{Float64}}}
    Κ::Vector{ElasticMatrix{Float64,Vector{Float64}}}
    CellCount::Vector{Int64}
end

#@kwdef struct CellParams_t
#    kₛ::Float64 = 0.0
#    η::Float64 = 1.0
#    kf::Float64
#    D::Float64 = 0.0
#    a::Float64 = 1.0
#end
