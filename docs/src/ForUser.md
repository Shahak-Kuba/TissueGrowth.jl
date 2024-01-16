```@meta
CurrentModule = TissueGrowth
```

# User Set Parameters
`N`: Number of cells in the simulation. [Integer]

`m`: Number of springs per cell. [Integer]

`R₀`: Radius or characteristic length of the initial shape. [Float]

`D`: Array of diffuision coefficients used to calculate cell stiffness. [Float]

`l₀`: Resting length of the spring per cell. [Float]

`kf`: Tissue production rate per cell. [Float]

`η`: Viscosity or damping coefficient per cell. [Float]

`growth_dir`: Direction of tissue growth (Options: 'inward' or 'outward'). [String]

`Tmax`: Total simulation time (in days). [Float]

`δt`: Time step for the numerical integration. [Float]

`btypes`: Types of boundary conditions (Options: "Sinewave" (1D only), "circle", "triangle", "square", "hex", "star", "cross"). [String]

`dist_type`: Distribution type for node placement ("Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"). [String]

`prolif`, `death`, `embed`: Boolean flags indicating cell behaviors. [Boolean]

`α`, `β`, `γ`: Parameters for cell behaviors. [Float]

`event_δt`: Time interval for periodic callback events. [Float]

The stochastic cell behaviour has not been implemented in the 1D simulation code therefore a set of parameters would not be set/ required to run those simulations.