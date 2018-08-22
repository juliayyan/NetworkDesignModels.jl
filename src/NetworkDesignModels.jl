module NetworkDesignModels

    import JuMP, Gurobi, TransitNetworks, LightGraphs, MathProgBase
    const TN = TransitNetworks

    include("coveragemodel.jl")
    include("utils.jl")
    include("masterproblem.jl")
    include("subproblem.jl")
    include("subproblem-model.jl")
    include("subproblem-utils.jl")
end
