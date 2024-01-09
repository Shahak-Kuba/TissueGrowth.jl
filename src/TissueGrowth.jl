module TissueGrowth

    include("MechanicalEqns.jl")
    include("CellBehaviours.jl")
    include("PoreBoundaries.jl")
    include("PlottingFncs1D.jl")
    include("PlottingFncs2D.jl")
    include("Simulation_1D.jl")
    include("Simulation_2D.jl")
    include("TissueGrowthODEproblem.jl")
    include("GeometrySolvers.jl")
    include("Misc.jl")
    include("DataStructs.jl")
    include("PostSimulation.jl")
    include("ModifierFncs.jl")

end
