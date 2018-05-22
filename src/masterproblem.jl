mutable struct MasterProblem
    np::TN.TransitNetworkProblem
    linelist::Vector{Vector{Int}}
    commutelines::Dict{Tuple{Int,Int},Vector{Int}}
    costs::Vector{Float64}
    model::JuMP.Model
    budget::JuMP.Variable
    x::Vector{JuMP.Variable}
    θ::JuMP.JuMPDict{JuMP.Variable}
    choseline::JuMP.JuMPDict{JuMP.ConstraintRef}
    bcon::JuMP.ConstraintRef
    xub::Array{JuMP.ConstraintRef}
    choseub::JuMP.JuMPDict{JuMP.ConstraintRef}
    solver
end

nonzerodests(np::TN.TransitNetworkProblem, u::Int) =
    find(np.odmatrix[u,:] .> 0)

linecost(np::TN.TransitNetworkProblem, line::Vector{Int}) = 
    sum(TN.haversinedistance(np, line[i], line[i+1])
        for i in 1:length(line)-1)

function addline!(
    np::TN.TransitNetworkProblem,
    commutelines::Dict{Tuple{Int,Int},Vector{Int}},
    line::Vector{Int},
    lineindex::Int)
    for u in line, v in line
        u == v && continue
        !in(v, nonzerodests(np,u)) && continue
        push!(commutelines[u,v], lineindex)
    end 
end 

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
        commutelines[u,v] = Int[]
    end 
    for l in 1:nlines 
        addline!(np, commutelines, linelist[l], l)
    end 

    # cost computation
    costs = [linecost(np, line) for line in linelist]

    rmp, budget, x, θ, choseline, bcon, xub, choseub = 
        mastermodel(np, linelist, commutelines, costs, 
                    solver)
    JuMP.fix(budget, initialbudget)
    
    MasterProblem(np, linelist, commutelines, costs,
        rmp, budget, x, θ, choseline, bcon, xub, choseub,
        solver)
end 

function mastermodel(
    np::TN.TransitNetworkProblem,
    linelist::Vector{Vector{Int}},
    commutelines::Dict{Tuple{Int,Int},Vector{Int}},
    costs::Vector{Float64},
    solver)

    const nlines = length(linelist)

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
        θ[u,v] <= sum(x[l] for l in commutelines[u,v])
    )

    # budget constraint
    JuMP.@constraint(rmp, bcon, dot(costs, x) <= budget)

    rmp, budget, x, θ, choseline, bcon, xub, choseub
end 

function addcolumn!(rmp::MasterProblem,
    line::Vector{Int})
    # recompute line information
    const nlinesold = length(rmp.linelist)
    push!(rmp.linelist, line)
    push!(rmp.linelist, reverse(line))
    addline!(rmp.np, rmp.commutelines, line, nlinesold+1)
    addline!(rmp.np, rmp.commutelines, reverse(line), nlinesold+2)
    push!(rmp.costs, linecost(rmp.np, line))
    push!(rmp.costs, linecost(rmp.np, reverse(line)))
    initialbudget = JuMP.getvalue(rmp.budget)

    rmp.model, rmp.budget, 
    rmp.x, rmp.θ, rmp.choseline, 
    rmp.bcon, rmp.xub, rmp.choseub = 
        mastermodel(rmp.np, 
            rmp.linelist, rmp.commutelines, rmp.costs, 
            rmp.solver)
    JuMP.fix(rmp.budget, initialbudget)
end 

function optimize(mp::MasterProblem, budget::Int)
    JuMP.fix(mp.budget, budget)
    JuMP.solve(mp.model)
    JuMP.getvalue(mp.x)
end