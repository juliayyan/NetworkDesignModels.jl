mutable struct DynamicNetworkCoverageModel
    np::TN.TransitNetworkProblem
    model::JuMP.Model
    budget::JuMP.Variable
    δ::JuMP.Variable
    x::Matrix{JuMP.Variable}
    θ::JuMP.JuMPDict{JuMP.Variable}
end

function DynamicNetworkCoverageModel(
        np::TN.TransitNetworkProblem,
        odmatrix::Array{<: Real, 3};
        initialbudget::Real = 0.0,
        solver = Gurobi.GurobiSolver(),
        budgettype::Symbol = :count,
        delta::Int = 0
    )
    nlines = length(np.lines)
    nperiods = size(odmatrix, 3)
    model = JuMP.Model(solver=solver)
    JuMP.@variable(model, x[l=1:nlines, t=1:nperiods], Bin)
    JuMP.@variable(model, budget >= 0)
    JuMP.@variable(model, δ >= 0)
    JuMP.@variable(model, θ[
            u=1:np.nstations, v=TN.dests(np,u), r=1:TN.nroutes(np,u,v), t=1:nperiods
        ] >= 0
    )
    JuMP.@objective(model, Max, sum(
        sum(odmatrix[u,v,t]*θ[u,v,r,t] for r in 1:TN.nroutes(np,u,v), t in 1:nperiods)
        for u in 1:np.nstations, v in TN.dests(np, u)
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
    JuMP.@constraint(model,
        sum(costs[l]*x[l,t] for l in 1:nlines, t in 1:nperiods) <= nperiods*budget
    )
    JuMP.@variable(model, z[l=1:nlines, t=1:(nperiods-1)] >= 0)
    JuMP.@constraint(model, [l=1:nlines, t=1:(nperiods-1)], z[l,t] >= x[l,t+1] - x[l,t])
    JuMP.@constraint(model, [l=1:nlines, t=1:(nperiods-1)], z[l,t] >= x[l,t] - x[l,t+1])
    JuMP.@constraint(model, [t=1:(nperiods-1)], sum(z[l,t] for l=1:nlines) <= δ)

    JuMP.fix(δ, delta)
    JuMP.fix(budget, initialbudget)

    DynamicNetworkCoverageModel(np, model, budget, δ, x, θ)
end

function optimize(nm::DynamicNetworkCoverageModel, budget::Real)
    JuMP.fix(nm.budget, budget)
    JuMP.solve(nm.model)
    round.(Int, JuMP.getvalue.(nm.x))
end

function optimize(nm::DynamicNetworkCoverageModel, budget::Real, delta::Int)
    JuMP.fix(nm.budget, budget); JuMP.fix(nm.δ, delta)
    JuMP.solve(nm.model)
    round.(Int, JuMP.getvalue.(nm.x))
end
