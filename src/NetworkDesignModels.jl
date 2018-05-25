module NetworkDesignModels

    import JuMP, Gurobi, TransitNetworks, LightGraphs
    const TN = TransitNetworks

    include("coveragemodel.jl")
    include("masterproblem.jl")
    include("subproblem.jl")

end
