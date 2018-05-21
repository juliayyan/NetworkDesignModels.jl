mutable struct MasterProblem
    np::TN.TransitNetworkProblem
    linelist::Vector{Vector{Int}}
    model::JuMP.Model
    budget::JuMP.Variable
    x::Vector{JuMP.Variable}
    θ::JuMP.JuMPDict{JuMP.Variable}
    choseline::JuMP.JuMPDict{JuMP.ConstraintRef}
    bcon::JuMP.ConstraintRef
    xub::Array{JuMP.ConstraintRef}
    choseub::JuMP.JuMPDict{JuMP.ConstraintRef}
end

nonzerodests(np::TN.TransitNetworkProblem, u::Int) =
    find(np.odmatrix[u,:] .> 0)

function MasterProblem(
    np::TN.TransitNetworkProblem;
    initialbudget::Int = 0,
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    linelist::Vector{Vector{Int}} = np.lines
    )
    
    const nlines = length(linelist)

    # (u,v) --> lines that connect u and v
    commutelines = Dict{Tuple{Int,Int},Vector{Int}}()
    for u in 1:np.nstations, v in nonzerodests(np,u)
        commutelines[(u,v)] = Int[]
    end 
    for l in 1:nlines 
        for u in linelist[l], v in linelist[l]
            u == v && continue
            !in(v, nonzerodests(np,u)) && continue
            push!(commutelines[(u,v)], l)
        end 
    end 

    rmp = JuMP.Model(solver=solver)
    JuMP.@variable(rmp, x[l=1:nlines] >= 0)
    JuMP.@variable(rmp, budget)
    JuMP.@variable(rmp,
        θ[u=1:np.nstations,v=nonzerodests(np,u)] >= 0
    )
    
    # maximize ridership
    JuMP.@objective(rmp, 
        Max, 
        sum(sum(np.odmatrix[u,v]*θ[u,v] for v in nonzerodests(np,u)) 
                                        for u in 1:np.nstations)
    )

    # x constraints
    JuMP.@constraint(rmp, 
        xub[l=1:nlines], 
        x[l] <= 1)
    for l in 1:nlines, l2 in (l+1):nlines
        if sort(linelist[l]) == sort(linelist[l2])
            JuMP.@constraint(rmp, x[l] == x[l2])
        end 
    end 

    # choice constraints
    JuMP.@constraint(rmp, 
        choseub[u=1:np.nstations, v=nonzerodests(np,u)],
        θ[u,v] <= 1
    )
    JuMP.@constraint(rmp,
        choseline[u=1:np.nstations, v=nonzerodests(np,u)],
        θ[u,v] <= sum(x[l] for l in commutelines[(u,v)])
    )

    # budget constraint
    costs = [sum(TN.haversinedistance(np, linelist[l][i], linelist[l][i+1])
             for i in 1:length(linelist[l])-1) 
             for l in 1:length(linelist)]
    JuMP.@constraint(rmp, bcon, dot(costs, x) <= budget)
    
    JuMP.fix(budget, initialbudget)
    
    MasterProblem(np, linelist, rmp, budget, x, θ, choseline, bcon, xub, choseub)
end 

function optimize(mp::MasterProblem, budget::Int)
    JuMP.fix(mp.budget, budget)
    JuMP.solve(mp.model)
    JuMP.getvalue(mp.x)
end