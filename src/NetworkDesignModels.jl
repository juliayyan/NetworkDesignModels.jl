module NetworkDesignModels

    import JuMP, Gurobi, TransitNetworks; const TN = TransitNetworks

    mutable struct NetworkCoverageModel
        np::TN.TransitNetworkProblem
        model::JuMP.Model
        budget::Vector{JuMP.Variable}
        x::Matrix{JuMP.Variable}
        θ::JuMP.JuMPDict{JuMP.Variable}
    end

    function NetworkCoverageModel(
            np::TN.TransitNetworkProblem;
            initialbudget::Vector{Int} = Int[0],
            solver = Gurobi.GurobiSolver(),
            budgettype::Symbol = :count
        )
        nlines = length(np.lines)
        model = JuMP.Model(solver=solver)
        JuMP.@variable(model, x[l=1:nlines,1], Bin)
        JuMP.@variable(model, budget)
        JuMP.@variable(model,
            θ[u=1:np.nstations,v=TN.dests(np,u), r=1:TN.nroutes(np,u,v)] >= 0
        )
        JuMP.@objective(model, Max, sum(sum(
                np.odmatrix[u,v]*sum(θ[u,v,r] for r in 1:TN.nroutes(np,u,v))
                for v in TN.dests(np, u)
            )
            for u in 1:np.nstations
        ))
        JuMP.@constraint(model, [u=1:np.nstations, v=TN.dests(np,u)],
            sum(θ[u,v,r] for r in 1:TN.nroutes(np,u,v)) <= 1
        )
        JuMP.@constraint(model,
            [u=1:np.nstations, v=TN.dests(np,u),
             r=1:TN.nroutes(np,u,v), s=TN.stages(np,u,v,r)],
            θ[u,v,r] <= sum(x[l] for l in 1:nlines if np.linesegments[l,s])
        )

        if budgettype == :count 
            costs = ones(nlines)
        elseif budgettype == :distance         
            costs = [sum(TransitNetworks.haversinedistance(np,np.lines[l][i],np.lines[l][i+1])
                         for i in 1:length(np.lines[l])-1) 
                     for l in 1:nlines]
        else
            error("Incompatible budgettype") 
        end 
        JuMP.@constraint(model, sum(costs[l]*x[l,1] for l in 1:nlines) <= budget)
        
        JuMP.fix(budget, initialbudget)

        NetworkCoverageModel(np, model, JuMP.Variable[budget], x, θ)
    end

    function optimize(nm::NetworkCoverageModel, budget::Int)
        JuMP.fix(nm.budget[1], budget)
        JuMP.solve(nm.model)
        round.(Int, JuMP.getvalue.(nm.x))
    end

    include("dynamicproblem.jl")
end
