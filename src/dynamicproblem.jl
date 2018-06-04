
function NetworkCoverageModel{T}(
        np::TN.TransitNetworkProblem{Array{T,3}};
        initialbudget::Vector{Int} = zeros(Int, size(np.odmatrix,3)),
        solver = Gurobi.GurobiSolver(),
        budgettype::Symbol = :count
    ) where {T <: Real}
    nlines = length(np.lines)
    nperiods = length(initialbudget)
    model = JuMP.Model(solver=solver)
    JuMP.@variable(model, x[l=1:nlines, t=1:nperiods], Bin)
    JuMP.@variable(model, budget[t=1:nperiods])
    JuMP.@variable(model, θ[
            u=1:np.nstations, v=TN.dests(np,u), r=1:TN.nroutes(np,u,v), t=1:nperiods
        ] >= 0
    )
    JuMP.@objective(model, Max, sum(sum(
            np.odmatrix[u,v]*sum(
                θ[u,v,r,t] for r in 1:TN.nroutes(np,u,v), t in 1:nperiods
            )
            for v in TN.dests(np, u)
        )
        for u in 1:np.nstations
    ))
    JuMP.@constraint(model, [u=1:np.nstations, v=TN.dests(np,u), t=1:nperiods],
        sum(θ[u,v,r,t] for r in 1:TN.nroutes(np,u,v)) <= 1
    )
    JuMP.@constraint(model,
        [u=1:np.nstations, v=TN.dests(np,u),
         r=1:TN.nroutes(np,u,v), s=TN.stages(np,u,v,r), t=1:nperiods],
        θ[u,v,r,t] <= sum(x[l,t] for l in 1:nlines if np.linesegments[l,s])
    )

    if budgettype == :count 
        costs = ones(nlines)
    elseif budgettype == :distance         
        costs = [sum(TransitNetworks.haversinedistance(np,np.lines[l][i],np.lines[l][i+1])
                     for i in 1:length(np.lines[l])-1) 
                 for l in 1:nlines]
    else
        error("Invalid budget type $budgettype. Use `:count` or `:distance`.") 
    end 
    JuMP.@constraint(model, [t=1:nperiods],
        sum(costs[l]*x[l,t] for l in 1:nlines) <= budget[t]
    )
    
    for t in 1:nperiods; JuMP.fix(budget[t], initialbudget[t]) end

    NetworkCoverageModel(np, model, budget, x, θ)
end

function optimize(nm::NetworkCoverageModel, initialbudget::Vector{Int})
    for t in 1:nperiods; JuMP.fix(nm.budget[t], initialbudget[t]) end
    JuMP.solve(nm.model)
    round.(Int, JuMP.getvalue.(nm.x))
end
