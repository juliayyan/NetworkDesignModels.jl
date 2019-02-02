module NetworkDesignModels

    import JuMP, Gurobi, TransitNetworks, LightGraphs, MathProgBase
    using LinearAlgebra
    TN = TransitNetworks

    include("coveragemodel.jl")
    include("utils.jl")
    include("masterproblem.jl")
    include("subproblem.jl")
    include("subproblem-traveltimes.jl")
    include("subproblem-model.jl")
    include("subproblem-utils.jl")
end
