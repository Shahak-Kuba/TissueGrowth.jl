var documenterSearchIndex = {"docs":
[{"location":"Example.html#Example-2D-simulation","page":"Example 2D simulation","title":"Example 2D simulation","text":"","category":"section"},{"location":"Example.html#Non-Stochastic-Simulation-Code:","page":"Example 2D simulation","title":"Non-Stochastic Simulation Code:","text":"","category":"section"},{"location":"Example.html","page":"Example 2D simulation","title":"Example 2D simulation","text":"using TissueGrowth\n\n# set random seed number for reproducability \nseed = 88\n\n# setting up simulation parameters\nN = 180 # number of cells\nm = 1 # number of springs per cell\nR₀ = 1.0  # shape radius\nD = [0.0001, 0.0075, 0.015]\nl₀ = 1.0\nkf = 0.0008\nη = 1.0 \ngrowth_dir = \"inward\" # Options: \"inward\", \"outward\"\nTmax = 21.0 # days\nδt = 0.01\nbtypes = [\"circle\"] #Options: [\"circle\", \"triangle\", \"square\", \"hex\", \"star\",\"cross\"]\ndist_type = \"Linear\" #Options: [\"Linear\", \"sigmoid\", \"2sigmoid\", \"exp\",  \"sine\", \"cosine\", \"quad\", \"cubic\"]\n\n# Cell Behaviours\nprolif = false;    death = false;  embed = false;\nα = 0.0001;        β = 0.001;      γ = 0.01;\nevent_δt = δt\n\n# 2D simulations \nsols2D = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,\n            prolif, death, embed, α, β, γ, event_δt, seed);\n\ncmap = :jet\n\ngeo = 1 # indexing through `btypes`\ndiffusivity = 1 # indexing through `D`\nTissueGrowth.plotResults2D(sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].Density, cmap)\n","category":"page"},{"location":"Example.html#Stochastic-Simulation-Code:","page":"Example 2D simulation","title":"Stochastic Simulation Code:","text":"","category":"section"},{"location":"Example.html","page":"Example 2D simulation","title":"Example 2D simulation","text":"using TissueGrowth\n\n# set random seed number for reproducability \nseed = 88\n\n# setting up simulation parameters\nN = 180 # number of cells\nm = 1 # number of springs per cell\nR₀ = 1.0  # shape radius\nD = [0.0001, 0.0075, 0.015]\nl₀ = 1.0\nkf = 0.0008\nη = 1.0 \ngrowth_dir = \"inward\" # Options: \"inward\", \"outward\"\nTmax = 21.0 # days\nδt = 0.01\nbtypes = [\"circle\"] #Options: [\"circle\", \"triangle\", \"square\", \"hex\", \"star\",\"cross\"]\ndist_type = \"Linear\" #Options: [\"Linear\", \"sigmoid\", \"2sigmoid\", \"exp\",  \"sine\", \"cosine\", \"quad\", \"cubic\"]\n\n# Cell Behaviours\nprolif = true;     death = true;    embed = false;\nα = 0.0001;        β = 0.001;      γ = 0.01;\nevent_δt = δt\n\n# 2D simulations \nsols2D = TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,\n            prolif, death, embed, α, β, γ, event_δt, seed);\n\ncmap = :jet\n\ngeo = 1 # indexing through `btypes`\ndiffusivity = 1 # indexing through `D`\nTissueGrowth.plotResults2D(sols2D[diffusivity][geo].u, sols2D[diffusivity][geo].Density, cmap)\n","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"CurrentModule = TissueGrowth","category":"page"},{"location":"ForUser.html#User-Set-Parameters","page":"User Set Parameters","title":"User Set Parameters","text":"","category":"section"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"N: Number of cells in the simulation. [Integer]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"m: Number of springs per cell. [Integer]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"R₀: Radius or characteristic length of the initial shape. [Float]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"D: Array of diffuision coefficients used to calculate cell stiffness. [Float]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"l₀: Resting length of the spring per cell. [Float]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"kf: Tissue production rate per cell. [Float]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"η: Viscosity or damping coefficient per cell. [Float]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"growth_dir: Direction of tissue growth (Options: 'inward' or 'outward'). [String]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"Tmax: Total simulation time (in days). [Float]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"δt: Time step for the numerical integration. [Float]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"btypes: Types of boundary conditions (Options: \"Sinewave\" (1D only), \"circle\", \"triangle\", \"square\", \"hex\", \"star\", \"cross\"). [String]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"dist_type: Distribution type for node placement (\"Linear\", \"sigmoid\", \"2sigmoid\", \"exp\",  \"sine\", \"cosine\", \"quad\", \"cubic\"). [String]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"prolif, death, embed: Boolean flags indicating cell behaviors. [Boolean]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"α, β, γ: Parameters for cell behaviors. [Float]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"event_δt: Time interval for periodic callback events. [Float]","category":"page"},{"location":"ForUser.html","page":"User Set Parameters","title":"User Set Parameters","text":"The stochastic cell behaviour has not been implemented in the 1D simulation code therefore a set of parameters would not be set/ required to run those simulations.","category":"page"},{"location":"Equations.html#Simulation-functions-and-code-description","page":"Simulation functions and code description","title":"Simulation functions and code description","text":"","category":"section"},{"location":"Equations.html#Index","page":"Simulation functions and code description","title":"Index","text":"","category":"section"},{"location":"Equations.html","page":"Simulation functions and code description","title":"Simulation functions and code description","text":"","category":"page"},{"location":"Equations.html#Data-Structure","page":"Simulation functions and code description","title":"Data Structure","text":"","category":"section"},{"location":"Equations.html","page":"Simulation functions and code description","title":"Simulation functions and code description","text":"TissueGrowth.SimResults_t","category":"page"},{"location":"Equations.html#TissueGrowth.SimResults_t","page":"Simulation functions and code description","title":"TissueGrowth.SimResults_t","text":"SimResults_t\n\nA data structure representing the results of a simulation which can then be used for visualisation with the plotting code PlottingFncs1D.jl and PlottingFncs2D.jl.\n\nThis structure contains various fields that store the results and relevant data from a simulation process. It includes the type of simulation, time steps, state vectors, forces, densities, and other relevant quantities.\n\nFields\n\nbtype: A String indicating the type of simulation or boundary condition.\nt: A Vector{Float64} representing the time steps at which the simulation results are recorded.\nu: A Vector of ElasticMatrix{Float64,Vector{Float64}} representing the state of the system at each time step.\n∑F: A Vector of ElasticVector{Float64,Vector{Float64}}, representing the sum of forces at each time step.\nDensity: A Vector of ElasticMatrix{Float64,Vector{Float64}} representing the density of the system at each time step.\nVₙ: A Vector{Vector{Float64}} representing the normal velocity at each time step.\nΩ: A Vector{Float64} representing the void area within the pore at each time step.\nψ: A Vector of ElasticMatrix{Float64,Vector{Float64}}, representing the stress of the springs/ cells at each time step.\nΚ: A Vector of ElasticMatrix{Float64,Vector{Float64}}, representing the approximate curvature of the moving boundary at each time step.\nCellCount: A Vector{Int64} representing the amount of active cells at a give time t.\n\n\n\n\n\n","category":"type"},{"location":"Equations.html#1D-simulation-code:","page":"Simulation functions and code description","title":"1D simulation code:","text":"","category":"section"},{"location":"Equations.html","page":"Simulation functions and code description","title":"Simulation functions and code description","text":"TissueGrowth.sim1D()","category":"page"},{"location":"Equations.html#TissueGrowth.sim1D-Tuple{}","page":"Simulation functions and code description","title":"TissueGrowth.sim1D","text":"sim1D()\n\nExecute a series of 1D mechanical relaxation simulations.\n\nThis function sets up and runs a series of 1D simulations for different stiffness values. It iterates over an array of stiffness coefficients, sets up the corresponding ODE problem for each case, solves it, and collects the results.\n\nSimulation Parameters\n\nN: Number of cells in the simulation.\nm: Number of springs per cell.\nM: Total number of springs along the interface.\nR₀: Radius or characteristic length of the initial shape.\nD: Array of diffuision coefficients used to calculate cell stiffness.\nl₀: Resting length of the spring per cell.\nkf: Tissue production rate per cell.\nη: Viscosity or damping coefficient per cell.\nTmax: Total simulation time (in days).\nδt: Time step for the Euler integration.\ngrowth_dir: Growth direction for the simulation ('1D').\nbtype: Type of boundary condition ('SineWave').\nsavetimes: Time points at which to save the simulation results.\n\nReturns\n\nA vector of SimResults_t objects, each representing the simulation results for a different stiffness coefficient.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#2D-simulation-code:","page":"Simulation functions and code description","title":"2D simulation code:","text":"","category":"section"},{"location":"Equations.html","page":"Simulation functions and code description","title":"Simulation functions and code description","text":"TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,prolif,death,embed,α,β,γ,event_δt,seed)","category":"page"},{"location":"Equations.html#TissueGrowth.sim2D-NTuple{20, Any}","page":"Simulation functions and code description","title":"TissueGrowth.sim2D","text":"sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,prolif,death,embed,α,β,γ,event_δt,seed)\n\nExecute a series of 2D mechanical relaxation simulations.\n\nThis function sets up and runs 2D simulations for different stiffness coefficients and boundary types. It iterates over arrays of stiffness coefficients and boundary types, sets up the corresponding ODE problem for each case, solves it using a specified numerical method with periodic callbacks, and collects the results.\n\nSimulation Parameters\n\nN: Number of cells in the simulation.\nm: Number of springs per cell.\nR₀: Radius or characteristic length of the initial shape.\nD: Array of diffuision coefficients used to calculate cell stiffness.\nl₀: Resting length of the spring per cell.\nkf: Tissue production rate per cell.\nη: Viscosity or damping coefficient per cell.\ngrowth_dir: Direction of tissue growth ('inward' or 'outward').\nTmax: Total simulation time (in days).\nδt: Time step for the numerical integration.\nbtypes: Types of boundary conditions (e.g., 'circle', 'triangle').\ndist_type: Distribution type for node placement (e.g., 'Linear', 'sigmoid').\nprolif, death, embed: Boolean flags indicating cell behaviors.\nα, β, γ: Parameters for cell behaviors.\nevent_δt: Time interval for periodic callback events.\nseed: Seed number for reproducability with stochastic simulations.\n\nReturns\n\nA vector of vectors of SimResults_t objects. Each inner vector represents the simulation results for different boundary types under a specific stiffness coefficient.\n\nExample\n\nall_results = sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,\n                        prolif, death, embed, α, β, γ, event_δt, seed);\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#ODE-Problem:","page":"Simulation functions and code description","title":"ODE Problem:","text":"","category":"section"},{"location":"Equations.html","page":"Simulation functions and code description","title":"Simulation functions and code description","text":"TissueGrowth.ODE_fnc_1D!(du,u,p,t) \nTissueGrowth.ODE_fnc_2D!(du,u,p,t) ","category":"page"},{"location":"Equations.html#TissueGrowth.ODE_fnc_1D!-NTuple{4, Any}","page":"Simulation functions and code description","title":"TissueGrowth.ODE_fnc_1D!","text":"ODE_fnc_1D!(du, u, p, t)\n\nDefine the ODE system for 1D mechanical relaxation with initial conditions. This function computes the derivatives du based on the current state u and parameters p. In this system there is a periodic boundary condition such that x ∈ [0,2π] with i spring boundary nodes and u₋₁ = uₙ, uₙ₊₁ = u₁\n\nArguments\n\ndu: Array to store the derivatives of u.\nu: Current state array.\np: Parameters tuple (N, kₛ, η, kf, l₀, δt, growth_dir).\nt: Current time.\n\nDescription\n\nCalculates the mechanical relaxation and normal velocity in a 1D system with periodic boundary conditions. The derivatives are based on spring forces and mechanical properties defined in p.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.ODE_fnc_2D!-NTuple{4, Any}","page":"Simulation functions and code description","title":"TissueGrowth.ODE_fnc_2D!","text":"ODE_fnc_2D!(du, u, p, t)\n\nDefine the ODE system for 2D mechanical relaxation with initial conditions. This function computes the derivatives du based on the current state u and parameters p.\n\nArguments\n\ndu: Array to store the derivatives of u.\nu: Current position array of all spring boundary nodes.\np: Parameters tuple (N, kₛ, η, kf, l₀, δt, growth_dir).\nt: Current time.\n\nDescription\n\nCalculates the mechanical relaxation and normal velocity in a 1D system with periodic boundary conditions. The derivatives are based on spring forces and mechanical properties defined in p.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#Equations-for-ODE-problem:","page":"Simulation functions and code description","title":"Equations for ODE problem:","text":"","category":"section"},{"location":"Equations.html","page":"Simulation functions and code description","title":"Simulation functions and code description","text":"TissueGrowth.δ(rᵢ₊₁, rᵢ)\nTissueGrowth.ρ(rᵢ₊₁, rᵢ)\nTissueGrowth.τ(rᵢ₊₁, rᵢ₋₁)\nTissueGrowth.n(rᵢ₊₁, rᵢ₋₁, type)\nTissueGrowth.Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)\nTissueGrowth.Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)\nTissueGrowth.Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt, type)\nTissueGrowth.κ(rᵢ₋₁, rᵢ, rᵢ₊₁)","category":"page"},{"location":"Equations.html#TissueGrowth.δ-Tuple{Any, Any}","page":"Simulation functions and code description","title":"TissueGrowth.δ","text":"δ(rᵢ₊₁, rᵢ)\n\nCalculate the Euclidean distance between two points rᵢ₊₁ and rᵢ.\n\nArguments\n\nrᵢ₊₁: The first point in space.\nrᵢ: The second point in space.\n\nReturns\n\nThe Euclidean distance between the two points.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.ρ-Tuple{Any, Any}","page":"Simulation functions and code description","title":"TissueGrowth.ρ","text":"ρ(rᵢ₊₁, rᵢ)\n\nCalculate the reciprocal of the distance (interpreted as density) between two points rᵢ₊₁ and rᵢ.\n\nArguments\n\nrᵢ₊₁: The first point in space.\nrᵢ: The second point in space.\n\nReturns\n\nThe reciprocal of the distance between the two points.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.τ-Tuple{Any, Any}","page":"Simulation functions and code description","title":"TissueGrowth.τ","text":"τ(rᵢ₊₁, rᵢ₋₁)\n\nCalculate the unit tangent vector between two neighboring points rᵢ₊₁ and rᵢ₋₁.\n\nArguments\n\nrᵢ₊₁: The point after the central point in space.\nrᵢ₋₁: The point before the central point in space.\n\nReturns\n\nThe unit tangent vector between the two points.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.n-Tuple{Any, Any, Any}","page":"Simulation functions and code description","title":"TissueGrowth.n","text":"n(rᵢ₊₁, rᵢ₋₁, type)\n\nCalculate the unit normal vector at a point rᵢ between two neighboring points rᵢ₊₁ and rᵢ₋₁. The orientation of the normal vector depends on the specified type.\n\nArguments\n\nrᵢ₊₁: The point after the central point in space.\nrᵢ₋₁: The point before the central point in space.\ntype: A string specifying the orientation of the normal vector, either \"inward\" or any other value for outward orientation.\n\nReturns\n\nThe unit normal vector at the point.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.Fₛ⁺-NTuple{5, Any}","page":"Simulation functions and code description","title":"TissueGrowth.Fₛ⁺","text":"Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)\n\nCalculate the spring force (Nonlinear) for mechanical relaxation in the positive direction.\n\nArguments\n\nrᵢ: The current point in space.\nrᵢ₊₁: The point after the current point in space.\nrᵢ₋₁: The point before the current point in space.\nkₛ: Spring coefficient.\nl₀: Resting length of the spring.\n\nReturns\n\nThe spring force in the positive direction.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.Fₛ⁻-NTuple{5, Any}","page":"Simulation functions and code description","title":"TissueGrowth.Fₛ⁻","text":"Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)\n\nCalculate the spring force (Nonlinear) for mechanical relaxation in the negative direction.\n\nArguments\n\nrᵢ: The current point in space.\nrᵢ₊₁: The point after the current point in space.\nrᵢ₋₁: The point before the current point in space.\nkₛ: Spring coefficient.\nl₀: Resting length of the spring.\n\nReturns\n\nThe spring force in the negative direction.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.Vₙ-NTuple{6, Any}","page":"Simulation functions and code description","title":"TissueGrowth.Vₙ","text":"Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt, type)\n\nCalculate the normal velocity of the interface such that Vₙ is proportional to ρ. The direction of the normal vector is determined by type.\n\nArguments\n\nrᵢ₋₁: The point before the current point in space.\nrᵢ: The current point in space.\nrᵢ₊₁: The point after the current point in space.\nkf: The amount of tissue produced per unit area per unit time.\nδt: The time step.\ntype: A string specifying the orientation of the normal vector, either \"inward\" or any other value for outward orientation.\n\nReturns\n\nThe normal velocity of the interface.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.κ-Tuple{Any, Any, Any}","page":"Simulation functions and code description","title":"TissueGrowth.κ","text":"κ(rᵢ₋₁, rᵢ, rᵢ₊₁)\n\nApproximate the curvature of a shape using the Menger method. This method is based on the areas of triangles formed by consecutive triplets of points.\n\nArguments\n\nrᵢ₋₁: The point before the current point in space.\nrᵢ: The current point in space.\nrᵢ₊₁: The point after the current point in space.\n\nReturns\n\nThe approximated curvature at the point rᵢ.\n\nReference\n\nAnoshkina, Elena V., Alexander G. Belyaev, and Hans-Peter Seidel. \"Asymtotic Analysis of Three-Point Approximations of Vertex Normals and Curvatures.\" VMV. 2002.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#Cell-behaviour-functions:","page":"Simulation functions and code description","title":"Cell behaviour functions:","text":"","category":"section"},{"location":"Equations.html","page":"Simulation functions and code description","title":"Simulation functions and code description","text":"TissueGrowth.P(event,ρ,α)\nTissueGrowth.A(event,ρ,β)\nTissueGrowth.E(event,ρ,γ)\nTissueGrowth.affect!(integrator)\nTissueGrowth.store_embed_cell_pos(pos)","category":"page"},{"location":"Equations.html#TissueGrowth.P-Tuple{Any, Any, Any}","page":"Simulation functions and code description","title":"TissueGrowth.P","text":"P(event, ρ, α)\n\nCalculate the probability of proliferation occurring, given the density ρ and the parameter α.\n\nIf the proliferation event is considered to occur (event is true), the probability is calculated as the product of ρ and α. Otherwise, a vector of zeros is returned.\n\nArguments\n\nevent: A boolean indicating whether the proliferation event is considered to occur.\nρ: A vector representing densities.\nα: Parameter for scaling the proliferation probability.\n\nReturns\n\nA vector representing the calculated probabilities for proliferation.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.A-Tuple{Any, Any, Any}","page":"Simulation functions and code description","title":"TissueGrowth.A","text":"A(event, ρ, β)\n\nCalculate the probability of apoptosis (cell death) occurring, given the density ρ and the parameter β.\n\nIf the apoptosis event is considered to occur (event is true), the probability is calculated as the product of ρ and β. Otherwise, a vector of zeros is returned.\n\nArguments\n\nevent: A boolean indicating whether the apoptosis event is considered to occur.\nρ: A vector representing densities.\nβ: Parameter for scaling the apoptosis probability.\n\nReturns\n\nA vector representing the calculated probabilities for apoptosis.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.E-Tuple{Any, Any, Any}","page":"Simulation functions and code description","title":"TissueGrowth.E","text":"E(event, ρ, γ)\n\nCalculate the probability of embedding occurring, given the density ρ and the parameter γ.\n\nIf the embedding event is considered to occur (event is true), the probability is calculated as the product of ρ and γ. Otherwise, a vector of zeros is returned.\n\nArguments\n\nevent: A boolean indicating whether the embedding event is considered to occur.\nρ: A vector representing densities.\nγ: Parameter for scaling the embedding probability.\n\nReturns\n\nA vector representing the calculated probabilities for embedding.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.affect!-Tuple{Any}","page":"Simulation functions and code description","title":"TissueGrowth.affect!","text":"affect!(integrator)\n\nUpdate the state of integrator based on probabilistic cellular events.\n\nThis function modifies integrator in place. It uses the parameters and state from integrator to compute probabilities for different cellular events: proliferation (prolif), death (death), and embedding (embed). Based on these probabilities and random draws, it updates the state of integrator.\n\nArguments\n\nintegrator: The integrator object containing the current state and parameters. The parameters are expected to be a tuple containing:\n\nReturns\n\nnothing. The function modifies integrator in place.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#TissueGrowth.store_embed_cell_pos-Tuple{Any}","page":"Simulation functions and code description","title":"TissueGrowth.store_embed_cell_pos","text":"store_embed_cell_pos(pos)\n\nStores the position of embedded cells into 🥔 array.\n\nThis function inserts the position vector pos of an embedded cell into the 🥔 array.\n\nArguments\n\npos: position vector of the cell that is being embedded into the tissue should contain m+1 values given m springs in the simulation\n\nReturns\n\nnothing. The function modifies 🥔 in place.\n\n\n\n\n\n","category":"method"},{"location":"Equations.html#Custom-modifier-functions","page":"Simulation functions and code description","title":"Custom modifier functions","text":"","category":"section"},{"location":"Equations.html","page":"Simulation functions and code description","title":"Simulation functions and code description","text":"Base.insert!\nBase.deleteat!","category":"page"},{"location":"Equations.html#Base.insert!","page":"Simulation functions and code description","title":"Base.insert!","text":"Base.insert!(A::ElasticMatrix, col, xy)\n\nInsert a two-element vector xy into a specific column col of an ElasticMatrix A.\n\nThis function modifies the ElasticMatrix in place. It inserts the elements of xy into the underlying data of A at the positions corresponding to the specified column.\n\nArguments\n\nA::ElasticMatrix: The ElasticMatrix object to be modified.\ncol: The column index where the elements will be inserted.\nxy: A two-element vector containing the values to be inserted into the matrix.\n\nReturns\n\nThe modified ElasticMatrix A.\n\n\n\n\n\n","category":"function"},{"location":"Equations.html#Base.deleteat!","page":"Simulation functions and code description","title":"Base.deleteat!","text":"Base.deleteat!(A::ElasticMatrix, col)\n\nDelete elements from a specific column col of an ElasticMatrix A.\n\nThis function modifies the ElasticMatrix in place. It removes the elements of the specified column from the underlying data of A.\n\nArguments\n\nA::ElasticMatrix: The ElasticMatrix object to be modified.\ncol: The column index from which the elements will be removed.\n\nReturns\n\nThe modified ElasticMatrix A.\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"CurrentModule = TissueGrowth","category":"page"},{"location":"index.html#TissueGrowth.jl","page":"Home","title":"TissueGrowth.jl","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"TissueGrowth.jl is a open source project package developed to simulate the evoluition of biological tissue interface during tissue growth. Our model considers mechanical interactions between neighbouring cells and a secretion rate of new tissue material that is proportional to cell density. In this project we include cell proliferation, apoptosis (death) and embedment as stochastic processes. This model is solved using the DifferentialEquations.jl and ElasticArrays.jl packages. The solver uses RK4() or Euler() solving algorithms with a constant timestep Delta t. The stochastic cell behaviour is implemented using a PeriodicCallback and occurs every delta t_textevent period of time. This package offers 1D and 2D simulations. In our case 1D simulations imply the evolution of a line segment which has periodic boundaries and 2D is a closed domain. In the case of 2D simulations you can choose an inward or outward growth simulation by setting the appropriate parameters.","category":"page"},{"location":"index.html#Model-Description","page":"Home","title":"Model Description","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"This cell-based mathematical model ","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Calculates the mechanical relaxation and normal velocity in a system with periodic boundary conditions. The derivatives are based on the discrete equation given by:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"fractextdmathbfuᵢtextdt = frac1ηbigg((mathbfFₛ - mathbfFₛ)cdotbmτbigg)bmτ + Vₙmathbfn ","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"where,","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"mathbfFₛ = f(lVert mathbfuᵢ₁ - mathbfuᵢ rVert) fracmathbfuᵢ₁ - mathbfuᵢlVert mathbfuᵢ₁ - mathbfuᵢ rVert hspace05cm mathbfFₛ = f(lVert mathbfuᵢ - mathbfuᵢ₁ rVert) fracmathbfuᵢ - mathbfuᵢ₁lVert mathbfuᵢ - mathbfuᵢ₁ rVert","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"given f(lVert mathbfuᵢ₁ - mathbfuᵢ rVert) is the restoration force function and Vₙ is the velocity in the normal direction.","category":"page"},{"location":"index.html#Current-models-of-Mechanical-Relaxation-and-Normal-Velocity","page":"Home","title":"Current models of Mechanical Relaxation and Normal Velocity","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"In this current version of the simulation code we use a nonlinear resorting force given by,","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"f(lVert mathbfuᵢ₁ - mathbfuᵢ rVert) = kₛ bigg(frac1l_0 - frac1lVert mathbfuᵢ₁ - mathbfuᵢ rVert bigg)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"and a normal velocity which is proportional to density such that,","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Vₙ = k_frho_i hspace05cm rho_i = frac1lVert mathbfuᵢ₁ - mathbfuᵢ rVert","category":"page"}]
}
