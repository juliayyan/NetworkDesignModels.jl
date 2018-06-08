mutable struct MasterProblem
    np::TN.TransitNetworkProblem
    linelist::Vector{Vector{Int}}
    commutelines::Vector{Dict{Tuple{Int,Int},Any}}
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

function MasterProblem(
    np::TN.TransitNetworkProblem;
    initialbudget::Int = 0,
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    linelist::Vector{Vector{Int}} = uniquelines(np.lines),
    nlegs::Int = 1
    )
    
    const nlines = length(linelist)
    linelistcopy = Vector{Int}[]

    # (u,v) --> lines that connect u and v
    @assert nlegs <= 2
    commutelines = fill(Dict{Tuple{Int,Int},Any}(), nlegs)
    for u in 1:np.nstations, v in nonzerodests(np,u)
        commutelines[1][u,v] = Int[]
        if nlegs == 2 
            commutelines[2][u,v] = Tuple{Int,Int}[]
        end
    end 
    for l in 1:nlines 
        addline!(np, linelistcopy, commutelines, linelist[l])
    end 

    # cost computation
    costs = [linecost(np, line) for line in linelistcopy]

    rmp, budget, x, θ, choseline, bcon, choseub = 
        mastermodel(np, linelistcopy, commutelines, costs, 
                    solver)
    JuMP.fix(budget, initialbudget)
    
    MasterProblem(np, linelistcopy, commutelines, costs,
        rmp, budget, x, θ, choseline, bcon, choseub,
        solver)
end 

function mastermodel(
    np::TN.TransitNetworkProblem,
    linelist::Vector{Vector{Int}},
    commutelines::Vector{Dict{Tuple{Int,Int},Any}},
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
    addline!(rmp.np, rmp.linelist, rmp.commutelines, line)
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

# adds line to `oldlines` and `commutelines`
function addline!(
    np::TN.TransitNetworkProblem,
    oldlines::Vector{Vector{Int}},
    commutelines::Vector{Dict{Tuple{Int,Int},Any}},
    line::Vector{Int};
    maxdot = 0.5)
    l1 = length(oldlines)+1
    # single-leg commutes
    for u in line, v in line
        u == v && continue
        !in(v, nonzerodests(np,u)) && continue
        push!(commutelines[1][u,v], l1)
    end 
    # two-leg commutes
    if length(commutelines) == 2
        for l2 in 1:length(oldlines)
            line2 = oldlines[l2]
            xfrstns = intersect(line,line2)
            length(xfrstns) == 0 && continue
            stns1 = setdiff(line , xfrstns)
            stns2 = setdiff(line2, xfrstns)
            for u in stns1, v in stns2
                for w in xfrstns 
                    dir1 = dir(np,u,w)
                    dir2 = dir(np,w,v)
                    if dot(dir1,dir2)/norm(dir1)/norm(dir2) <= maxdot
                        haskey(commutelines[2], (u,v)) && 
                            push!(commutelines[2][(u,v)], (l1,l2))
                        haskey(commutelines[2], (v,u)) && 
                            push!(commutelines[2][(v,u)], (l2,l1))
                        break
                    end
                end        
            end
        end
    end
    push!(oldlines, line)
end 

function optimize(mp::MasterProblem, budget::Int)
    JuMP.fix(mp.budget, budget)
    JuMP.solve(mp.model)
    JuMP.getvalue(mp.x)
end
