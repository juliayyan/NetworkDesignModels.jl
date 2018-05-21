mutable struct NetworkCoverageModel
    np::TN.TransitNetworkProblem
    model::JuMP.Model
    budget::JuMP.Variable
    x::Vector{JuMP.Variable}
    θ::JuMP.JuMPDict{JuMP.Variable}
end

function NetworkCoverageModel(
        np::TN.TransitNetworkProblem;
        initialbudget::Int=0,
        solver = Gurobi.GurobiSolver(),
        budgettype::Symbol = :count
    )
    model = JuMP.Model(solver=solver)
    JuMP.@variable(model, x[l=1:length(np.lines)], Bin)
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
        θ[u,v,r] <= sum(x[l] for l in 1:length(np.lines) if np.linesegments[l,s])
    )

    if budgettype == :count 
        costs = ones(length(np.lines))
    elseif budgettype == :distance         
        costs = [sum(TransitNetworks.haversinedistance(np,np.lines[l][i],np.lines[l][i+1])
                     for i in 1:length(np.lines[l])-1) 
                 for l in 1:length(np.lines)]
    else
        error("Incompatible budgettype") 
    end 
    JuMP.@constraint(model, dot(costs, x) <= budget)
    
    JuMP.fix(budget, initialbudget)

    NetworkCoverageModel(np, model, budget, x, θ)
end

function optimize(nm::NetworkCoverageModel, budget::Int)
    JuMP.fix(nm.budget, budget)
    JuMP.solve(nm.model)
    round.(Int, JuMP.getvalue.(nm.x))
end