module NetworkDesignModels
    
    export DesignModel

    import DataFrames, Iterators, JuMP, Gurobi
    using TransitNetworks

    include("model/designmodel.jl")
    include("model/variables.jl")

end
