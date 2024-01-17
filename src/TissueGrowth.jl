module TissueGrowth

    include("Discrete/MechanicalEqns.jl")
    include("Discrete/CellBehaviours.jl")
    include("Discrete/PoreBoundariesV2.jl")
    include("Discrete/ProblemSetup.jl")
    include("Discrete/PlottingFncs1D.jl")
    include("Discrete/PlottingFncs2D.jl")
    include("Discrete/Simulation_1D.jl")
    include("Discrete/Simulation_2D.jl")
    include("Discrete/TissueGrowthODEproblem.jl")
    include("Discrete/GeometrySolvers.jl")
    include("Discrete/Misc.jl")
    include("Discrete/DataStructs.jl")
    include("Discrete/PostSimulation.jl")
    include("Discrete/ModifierFncs.jl")

end
