mutable struct MasterProblem
    # Network information
    np::TN.TransitNetworkProblem
    linelist::Vector{Vector{Int}}
    commutelines::Vector{Dict{Tuple{Int,Int},Any}}
    angleparam::Float64
    distparam::Float64
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
* `linelist`: A vector of the lines (assumed to be unique) in the transit
    network. Each line is a vector of ints corresponding to the stations on it.
* `nlegs`: The maximum number of lines being used for any given commute. So if
    `nlegs=2`, then only a single transfer is allowed.
* `angleparam`: -1 to 1. Lower values for sharper-angled transfers.
* `distparam`: >= 1.0.  Higher values allow longer commutes.
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
        angleparam::Float64 = -1.0, # deprecated by default
        distparam::Float64 = 1.5,
        solver = Gurobi.GurobiSolver(OutputFlag = 0),
        modeltype::Symbol = :lp
    )
    @assert in(modeltype, [:lp, :ip])
    @assert in(gridtype, [:latlong, :euclidean])
    @assert abs(angleparam) <= 1.0
    @assert distparam >= 1.0

    commutelines = allcommutelines(np, nlegs, linelist, angleparam, distparam, gridtype)
    costs = [linecost(np, line, gridtype) for line in linelist]
    rmp, budget, x, θ, choseline, bcon, choseub, pair1, pair2 = mastermodel(
        np, linelist, commutelines, costs, solver, modeltype
    )
    
    MasterProblem(
        np, linelist, commutelines, angleparam, distparam, costs, gridtype, rmp, budget,
        x, θ, choseline, bcon, choseub, pair1, pair2, solver, modeltype
    )
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
        JuMP.@constraint(rmp, pair1[p in pairs], aux[p] <= x[p[1]])
        JuMP.@constraint(rmp, pair2[p in pairs], aux[p] <= x[p[2]])
    else
        pair1 = pair2 = nothing
    end
    JuMP.@variable(rmp, budget)
    JuMP.@variable(rmp, θ[u=1:np.nstations,v=nonzerodests(np,u)])
    
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
            0 : sum(aux[pair] for pair in commutelines[2][u,v]))
    )
    # budget constraint
    JuMP.@constraint(rmp, bcon, dot(costs, x) <= budget)

    rmp, budget, x, θ, choseline, bcon, choseub, pair1, pair2
end 

"""
Adds a new column to the master problem.

It updates network information and creates a new model.
"""
function addcolumn!(
        rmp::MasterProblem,
        line::Vector{Int}
    )
    # recompute line information
    addline!(
        rmp.np, rmp.linelist, rmp.commutelines, line, 
        rmp.angleparam, rmp.distparam, rmp.gridtype
    )
    push!(rmp.costs, linecost(rmp.np, line,rmp.gridtype))
    initialbudget = JuMP.getvalue(rmp.budget)

    rmp.model, rmp.budget, rmp.x, rmp.θ, rmp.choseline, rmp.bcon, rmp.choseub,
        rmp.pair1, rmp.pair2 = mastermodel(
            rmp.np, rmp.linelist, rmp.commutelines, rmp.costs, rmp.solver,
            rmp.modeltype
        )
    JuMP.fix(rmp.budget, initialbudget)
end

"""
Returns a vector (with length `nlegs<=2`) of Dict{Tuple{Int,Int},Vector}, where

    commutelines[i][u,v] is a vector of the i lines that connect u to v.

In practice, we don't support i > 2, so 

    commutelines[1][u,v] is a Vector{Int}, and
    commutelines[2][u,v] is a Vector{Tuple{Int,Int}}.
"""
function allcommutelines(
        np::TN.TransitNetworkProblem,
        nlegs::Int,
        linelist::Vector{Vector{Int}},
        angleparam::Float64,
        distparam::Float64,
        gridtype::Symbol
    )
    @assert nlegs <= 2

    commutelines = [Dict{Tuple{Int,Int},Any}() for i in 1:nlegs]
    # 1. Initialize each (u,v) entry as an empty vector.
    for u in 1:np.nstations, v in nonzerodests(np,u)
        commutelines[1][u,v] = Int[]
        if nlegs == 2
            commutelines[2][u,v] = Tuple{Int,Int}[]
        end
    end
    # 2. Populate the (u,v) entries in commutelines.
    # 
    # We begin with an empty list, and iteratively add each line using the
    # addline!() method. At the end, we should recover the original linelist.
    linelistcopy = Vector{Int}[]
    for line in linelist
        linelistcopy = addline!(
            np, linelistcopy, commutelines, line, angleparam, distparam, gridtype
        )
    end
    @assert linelistcopy == linelist

    return commutelines
end

"""
Modify `oldlines` and `commutelines` in-place, adding information about `line`.
"""
function addline!(
        np::TN.TransitNetworkProblem,
        oldlines::Vector{Vector{Int}},
        commutelines::Vector{Dict{Tuple{Int,Int},Any}},
        line::Vector{Int},
        angleparam::Float64,
        distparam::Float64,
        gridtype::Symbol
    )
    l1 = length(oldlines) + 1 # We introduce a new `line` with index `l1`.
    # Single-leg commutes
    for u in line, v in line
        in(v, nonzerodests(np, u)) && push!(commutelines[1][u,v], l1)
    end
    # Two-leg commutes
    if length(commutelines) == 2
        for l2 in 1:length(oldlines)
            line2 = oldlines[l2]
            xfrstns = intersect(line, line2)
            length(xfrstns) == 0 && continue
            # u --line-> (w in xfrstns)
            #            (w in xfrstns) --line2-> v
            for u in setdiff(line, xfrstns), v in setdiff(line2, xfrstns)
                for w in xfrstns
                    if validtransfer(np, u, v, w, angleparam, distparam, gridtype)
                        pair = (min(l1, l2), max(l1, l2))
                        haskey(commutelines[2], (u, v)) &&
                            push!(commutelines[2][u,v], pair)
                        haskey(commutelines[2], (v, u)) &&
                            push!(commutelines[2][v,u], pair)
                        break
                    end
                end
            end
        end
    end
    # We update oldlines at the end to avoid conflicts with the earlier updates.
    push!(oldlines, line)
end

"""
Solves `mp::MasterProblem` with the given `budget`.

### Returns
The optimal solution `x` as a vector of floats.
"""
function optimize(mp::MasterProblem, budget::Real)
    JuMP.fix(mp.budget, budget)
    JuMP.solve(mp.model)
    JuMP.getvalue(mp.x)
end
