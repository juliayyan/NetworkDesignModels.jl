mutable struct MasterProblem
    # Network information
    np::TN.TransitNetworkProblem
    linelist::Vector{Vector{Int}}
    commutelines::Vector{Dict{Tuple{Int,Int},Any}}
    transferparam::Float64
    costs::Vector{Float64}
    gridtype::Symbol
    
    # Optimization information
    model::JuMP.Model
    budget::JuMP.Variable
    x::Vector{JuMP.Variable}
    θ::JuMP.JuMPDict{JuMP.Variable}
    choseline::JuMP.JuMPDict{JuMP.ConstraintRef}
    bcon::JuMP.ConstraintRef
    choseub::JuMP.JuMPDict{JuMP.ConstraintRef}
    pair1
    pair2
    solver
    modeltype::Symbol
end

"""
Returns a MasterProblem object.

Contains information on both the transit network and the optimization model.

### Keyword Arguments
* `gridtype`: `:latlong` or `:euclidean`.
* `linelist`:
* `nlegs`: The maximum number of lines being used for any given commute. So if
    `nlegs=2`, then only a single transfer is allowed.
* `transferparam`:
* `solver`: The solver being used to solve the problem.
* `modeltype`: `ip` or `lp`, to determine whether the relaxation is used.

### Quick Example
```
np = load("data/processed/networks/9b-transitnetwork.jld2", "keynetwork")
rmp = NetworkDesignModels.MasterProblem(np, nlegs = 1)
NetworkDesignModels.optimize(rmp, 5.0) # solves the problem with a budget of 5
```
"""
function MasterProblem(
        np::TN.TransitNetworkProblem;
        gridtype::Symbol = :latlong,
        linelist::Vector{Vector{Int}} = uniquelines(np.lines),
        nlegs::Int = 1,
        transferparam::Float64 = 0.0, # -1 to 1. lower allows sharper-angled transfers
        solver = Gurobi.GurobiSolver(OutputFlag = 0),
        modeltype::Symbol = :lp
    )
    @assert in(modeltype, [:lp, :ip])
    @assert in(gridtype, [:latlong, :euclidean])
    @assert nlegs <= 2
    
    # `commutelines` is a vector of Dict{Tuple{Int,Int},Vector}, where
    #
    #     commutelines[i][u,v] is a vector of the i lines that connect u to v.
    #
    # In practice, we don't support i > 2, so 
    #
    #     commutelines[1][u,v] is a Vector{Int}, and
    #     commutelines[2][u,v] is a Vector{Tuple{Int,Int}}.
    #
    commutelines = [Dict{Tuple{Int,Int},Any}() for i in 1:nlegs]
    # 1. We initialize each (u,v) entry in commutelines as an empty vector.
    for u in 1:np.nstations, v in nonzerodests(np,u)
        # (u,v) --> lines that connect u and v
        commutelines[1][u,v] = Int[]
        if nlegs == 2 
            commutelines[2][u,v] = Tuple{Int,Int}[]
        end
    end
    # 2. We populate the (u,v) entries in commutelines.
    linelistcopy = Vector{Int}[]
    for line in linelist
        addline!(np, linelistcopy, commutelines, line, transferparam,gridtype)
    end

    # cost computation
    costs = [linecost(np, line, gridtype) for line in linelistcopy]

    rmp, budget, x, θ, choseline, bcon, choseub, pair1, pair2 = 
        mastermodel(np, linelistcopy, commutelines, costs, 
                    solver, modeltype)
    
    MasterProblem(
        np, linelistcopy, commutelines, transferparam, costs, gridtype,
        rmp, budget, x, θ, choseline, bcon, choseub, pair1, pair2,
        solver, modeltype)
end 

"build base master problem model"
function mastermodel(
        np::TN.TransitNetworkProblem,
        linelist::Vector{Vector{Int}},
        commutelines::Vector{Dict{Tuple{Int,Int},Any}},
        costs::Vector{Float64},
        solver,
        modeltype::Symbol
    )

    const nlines = length(linelist)

    rmp = JuMP.Model(solver=solver)
    if modeltype == :lp
        JuMP.@variable(rmp, x[l=1:nlines] >= 0)
    elseif modeltype == :ip
        JuMP.@variable(rmp, x[l=1:nlines], Bin)
    end
    if length(commutelines) == 2
        pairs = unique(vcat(values(commutelines[2])...))
        JuMP.@variable(rmp, aux[pairs])
        JuMP.@constraint(rmp,
            pair1[p in pairs],
            aux[p] <= x[p[1]])
        JuMP.@constraint(rmp,
            pair2[p in pairs],
            aux[p] <= x[p[2]])
    else
        pair1 = pair2 = nothing
    end
    JuMP.@variable(rmp, budget)
    JuMP.@variable(rmp,
        θ[u=1:np.nstations,v=nonzerodests(np,u)]
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
        θ[u,v] <= sum(x[l] for l in commutelines[1][u,v]) +
        ((length(commutelines) == 1) || (length(commutelines[2][u,v]) == 0) ? 
            0 :
            sum(aux[pair] for pair in commutelines[2][u,v]))
    )
    
    # budget constraint
    JuMP.@constraint(rmp, bcon, dot(costs, x) <= budget)

    rmp, budget, x, θ, choseline, bcon, choseub, pair1, pair2
end 

"adds a new column to the master problem, updating network
 information and creating a new model"
function addcolumn!(
        rmp::MasterProblem,
        line::Vector{Int}
    )
    # recompute line information
    addline!(rmp.np, rmp.linelist, rmp.commutelines, line, 
        rmp.transferparam, rmp.gridtype)
    push!(rmp.costs, linecost(rmp.np, line,rmp.gridtype))
    initialbudget = JuMP.getvalue(rmp.budget)

    rmp.model, rmp.budget, 
    rmp.x, rmp.θ, rmp.choseline, 
    rmp.bcon, rmp.choseub, rmp.pair1, rmp.pair2 = 
        mastermodel(rmp.np, 
            rmp.linelist, rmp.commutelines, rmp.costs, 
            rmp.solver, rmp.modeltype)
    JuMP.fix(rmp.budget, initialbudget)
end 

"""
Modify `oldlines` and `commutelines` in-place, adding information about `line`.

### Remark
This is used as an utility function inside the constructor for MasterProblem,
and is not meant for ordinary usage.
"""
function addline!(
        np::TN.TransitNetworkProblem,
        oldlines::Vector{Vector{Int}},
        commutelines::Vector{Dict{Tuple{Int,Int},Any}},
        line::Vector{Int},
        transferparam::Float64,
        gridtype::Symbol
    )
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
                    if validtransfer(np,u,v,w,transferparam,gridtype)
                        pair = (min(l1,l2), max(l1,l2))
                        haskey(commutelines[2], (u,v)) && 
                            push!(commutelines[2][u,v], pair)
                        haskey(commutelines[2], (v,u)) && 
                            push!(commutelines[2][v,u], pair)
                        break
                    end
                end        
            end
        end
    end
    push!(oldlines, line)
end 

"""
Solves `mp::MasterProblem` with the given `budget`.

### Returns
The optimal solution `x` as a vector of floats.
"""
function optimize(mp::MasterProblem, budget::Float64)
    JuMP.fix(mp.budget, budget)
    JuMP.solve(mp.model)
    JuMP.getvalue(mp.x)
end
