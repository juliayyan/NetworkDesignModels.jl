mutable struct DesignModel
    np::TransitNetworkProblem
    model::JuMP.Model
    
    auxinfo::Dict{Symbol,Any}
    results::Dict{Symbol,Any}

    y::Vector{JuMP.Variable}
    ride::JuMP.JuMPDict{JuMP.Variable}
    open::JuMP.JuMPArray{JuMP.Variable,1,Tuple{Array{Array{Int64,1},1}}}

    function DesignModel(
        np::TransitNetworkProblem,
        model::JuMP.Model
    )
        new(np, model, Dict(), Dict())
    end
end

function DesignModel(
        np::TransitNetworkProblem,
        solver = Gurobi.GurobiSolver(OutputFlag=0),
        warmstart = zeros(0,0)
    )
    dm = DesignModel(np, JuMP.Model(solver=solver))

    # Variables
    dm.y    = designvariable(dm)
    dm.ride = ridershipvariable(dm)
    dm.open = openroutevariable(dm)

    # Constraints
    routeopenconstraint(dm)

    dm
end 
