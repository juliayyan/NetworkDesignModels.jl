module NetworkDesignModels

    import JuMP, Gurobi, TransitNetworks

    mutable struct NetworkCoverageModel
        np::TransitNetworks.TransitNetworkProblem
        model::JuMP.Model
        budget::JuMP.Variable
        x::Vector{JuMP.Variable}
        θ::JuMP.JuMPDict{JuMP.Variable}
    end

    function NetworkCoverageModel(
            np::TransitNetworks.TransitNetworkProblem;
            initialbudget::Int=0,
            solver = Gurobi.GurobiSolver(OutputFlag=0)
        )
        model = JuMP.Model(solver=solver)
        JuMP.@variable(model, x[l=1:length(np.lines)], Bin)
        JuMP.@variable(model, budget)
        JuMP.@variable(model,
            θ[u=1:np.nstations,v=TransitNetworks.dests(np,u),
              r=1:length(np.stage[u,v])] >= 0
        )
        JuMP.@objective(model, Max, sum(sum(
                np.odmatrix[u,v]*sum(θ[u,v,r] for r in 1:length(np.stage[u,v]))
                for v in TransitNetworks.dests(np, u)
            )
            for u in 1:np.nstations
        ))
        JuMP.@constraint(model, [u=1:np.nstations, v=TransitNetworks.dests(np,u)],
            sum(θ[u,v,r] for r in 1:length(np.stage[u,v])) <= 1
        )
        JuMP.@constraint(model,
            [u=1:np.nstations, v=TransitNetworks.dests(np,u),
             r=1:length(np.stage[u,v]), s=np.stage[u,v][r]],
            θ[u,v,r] <= sum(x[l] for l in 1:length(np.lines) if np.linesegments[l,s])
        )
        JuMP.@constraint(model, sum(x) <= budget)

        JuMP.fix(budget, initialbudget)

        NetworkCoverageModel(np, model, budget, x, θ)
    end

    function optimize(nm::NetworkCoverageModel, budget::Int)
        JuMP.fix(nm.budget, budget)
        JuMP.solve(nm.model)
        round.(Int, JuMP.getvalue.(nm.x))
    end

end
