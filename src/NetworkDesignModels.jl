module NetworkDesignModels

    import JuMP, Gurobi, TransitNetworks; const TN = TransitNetworks

    include("coveragemodel.jl")
    include("masterproblem.jl")
    include("subproblem.jl")

end
