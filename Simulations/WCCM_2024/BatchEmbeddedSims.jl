
using TissueGrowth
using Makie
using Printf
using BenchmarkTools

BatchSize = 1000

# See parameter approximation document
# Calculating kf
KF = 8784.2;
Tb = 28.46
l = 500;
Ω₀ = l^2
P = l*4
q₀ = 1/20; 
N = 50#Int(P*q₀) # number of cells
kf = KF/N
l_min = 5
l_max = 20


# setting up simulation parameters
m = 6 # number of springs per cell
R₀ = 75 #282.095  # shape radius μm
D = 0.00
kₛ = 1
Kₛ = 15
l₀ = 10.0
L₀ = ((l_max - l_min)/((kₛ/Kₛ)*((l_max^2 - l_min^2)/2 + l₀*(l_min - l_max)) - log(l_min/l_max)))
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 10 # days
δt = 0.01
btypes = ["PerturbedCircle"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m

## Cell Behaviours
prolif = false; death = false; embed = true;
α = 0.0;        β = 0.0;      Ot = 0.000625;
event_δt = δt

embedded_count_iteration_results = Vector{Int64}[]
Ω_iteration_results = Vector{Float64}[]
t = Vector{Float64}[]
all_solutions = Vector{TissueGrowth.SimResults_t}[]
all_embedded_cell_pos = Vector{Matrix{Float64}}[]


@time for iteration = 1:BatchSize
    seed = iteration
    sol, embedded_cell_pos, embedded_cell_count = TissueGrowth.GrowthSimulation(N,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 241);
    push!(all_solutions, sol)
    push!(all_embedded_cell_pos, embedded_cell_pos)
    push!(embedded_count_iteration_results, convert(Vector{Int64}, embedded_cell_count[1]))
    push!(Ω_iteration_results, sol[1].Ω[1] .- sol[1].Ω)
    if iteration == 1
        push!(t, sol[1].t)
    end
    println("Simulation $iteration / $BatchSize")
end

# converting vector of vectors into a Matrix
embedded_count_iteration_results_mat = reduce(vcat,embedded_count_iteration_results')
Ω_iteration_results_mat = reduce(vcat,Ω_iteration_results')
Ot_iteration_results_mat = embedded_count_iteration_results_mat ./ Ω_iteration_results_mat
Ot_iteration_results_mat[:,1] .= zeros(size(Ot_iteration_results_mat[:,1])) 
# calculating mins and max
min_Ot = minimum.(eachcol(Ot_iteration_results_mat))
max_Ot = maximum.(eachcol(Ot_iteration_results_mat))
σ_Ot = std.(eachcol(Ot_iteration_results_mat))

#Averaging Data
Ot_average = reduce(vcat,sum(Ot_iteration_results_mat,dims=1)./size(Ot_iteration_results_mat,1))
#Ω_average = sum(Ω_iteration_results)./size(Ω_iteration_results,1);

f = TissueGrowth.plotOtValueVsTime(t[1], Ot_average, Ot, min_Ot, max_Ot, m, σ_Ot)
save("WCCM_2024_PerturbedCircle_Batch_Plot_$BatchSize.png",f)
#f2 = TissueGrowth.plotOtValueVsTime(t[1], Ω_iteration_results[95], embedded_count_average, Ot/m)


geo = 1
diffusivity = 1

Density_cmap =  :cool #:rainbow1
Density_Range = (0.02,0.1)

#for iteration in 1:50
#    f = TissueGrowth.plotResults2D_embedded(all_solutions[iteration][1].u, all_solutions[iteration][1].Density, Density_cmap, Density_Range, "q [1/μm]", D, kf, (200,200), all_embedded_cell_pos[iteration], true)
#    save("WCCM_2024_Multi_embedded_$iteration.png",f)
#end
