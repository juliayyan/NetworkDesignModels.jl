mutable struct MasterProblem
    np::TN.TransitNetworkProblem
    linelist::Vector{Vector{Int}}
    commutelines::Vector{Dict{Tuple{Int,Int},Vector{Int}}}
    costs::Vector{Float64}
    model::JuMP.Model
    budget::JuMP.Variable
    x::Vector{JuMP.Variable}
    θ::JuMP.JuMPDict{JuMP.Variable}
    choseline::JuMP.JuMPDict{JuMP.ConstraintRef}
    bcon::JuMP.ConstraintRef
    choseub::JuMP.JuMPDict{JuMP.ConstraintRef}
    solver
end

function addline!(
    np::TN.TransitNetworkProblem,
    commutelines::Vector{Dict{Tuple{Int,Int},Vector{Int}}},
    line::Vector{Int},
    lineindex::Int)
    for u in line, v in line
        u == v && continue
        !in(v, nonzerodests(np,u)) && continue
        push!(commutelines[1][u,v], lineindex)
    end 
end 

function MasterProblem(
    np::TN.TransitNetworkProblem;
    initialbudget::Int = 0,
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    linelist::Vector{Vector{Int}} = uniquelines(np.lines),
    nlegs::Int = 1
    )
    
    const nlines = length(linelist)

    # (u,v) --> lines that connect u and v
    commutelines = fill(Dict{Tuple{Int,Int},Vector{Int}}(), nlegs)
    for u in 1:np.nstations, v in nonzerodests(np,u), i in 1:nlegs
        commutelines[i][u,v] = Int[]
    end 
    for l in 1:nlines 
        addline!(np, commutelines, linelist[l], l)
    end 

    # cost computation
    costs = [linecost(np, line) for line in linelist]

    rmp, budget, x, θ, choseline, bcon, choseub = 
        mastermodel(np, linelist, commutelines, costs, 
                    solver)
    JuMP.fix(budget, initialbudget)
    
    MasterProblem(np, linelist, commutelines, costs,
        rmp, budget, x, θ, choseline, bcon, choseub,
        solver)
end 

function mastermodel(
    np::TN.TransitNetworkProblem,
    linelist::Vector{Vector{Int}},
    commutelines::Vector{Dict{Tuple{Int,Int},Vector{Int}}},
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

    # choice constraints
    JuMP.@constraint(rmp, 
        choseub[u=1:np.nstations, v=nonzerodests(np,u)],
        θ[u,v] <= 1
    )
    JuMP.@constraint(rmp,
        choseline[u=1:np.nstations, v=nonzerodests(np,u)],
        θ[u,v] <= sum(x[l] for l in commutelines[1][u,v])
    )

    # budget constraint
    JuMP.@constraint(rmp, bcon, dot(costs, x) <= budget)

    rmp, budget, x, θ, choseline, bcon, choseub
end 

function addcolumn!(rmp::MasterProblem,
    line::Vector{Int})
    # recompute line information
    const nlinesold = length(rmp.linelist)
    push!(rmp.linelist, line)
    addline!(rmp.np, rmp.commutelines, line, nlinesold+1)
    push!(rmp.costs, linecost(rmp.np, line))
    initialbudget = JuMP.getvalue(rmp.budget)

    rmp.model, rmp.budget, 
    rmp.x, rmp.θ, rmp.choseline, 
    rmp.bcon, rmp.choseub = 
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
