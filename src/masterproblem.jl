"""
MasterOptions: A struct containing parameters for a MasterProblem

### Keyword Arguments
* `nlegs`: 1 (direct-route model) or 2 (single-transfer model)
* `nfreqs`: number of frequency levels to choose from
* `freqwts`: gamma in the paper; relative direct ridership coefficients
* `xfrwts`: lambda in the paper; relative transferring ridership coefficients
* `costwts`: rho in the paper; relative cost coefficients
* `constrainedg`: whether to constrain edges to be used at most once
* `ninitial`: the first `ninitial` lines in the `linelist` were used to start the master problem
* `angleparam`: old parameter; keep to -1.0
* `distparam`: Gamma in the paper; higher values allow longer commutes
"""
@with_kw mutable struct MasterOptions
    nlegs::Int                  = 1
    nfreqs::Int                 = 1
    freqwts::Vector{Float64}    = [1.0]
    xfrwts::Vector{Float64}     = [1.0]
    costwts::Vector{Float64}    = [1.0]
    constrainedg::Bool          = false
    ninitial::Int               = 0   
    angleparam::Float64         = -1.0 # unused by default
    distparam::Float64          = 1.5
    modeltype::Symbol           = :lp
end

"""
MasterProblem

* `np`: a TransitNetwork object that contains the network information
* `linelist`: a list of lines to choose to operate
* `commutelines`: has length `options.nlegs`.
                  The first entry maps a commute (u, v) to the lines that contain both stops
                  The second entry is unused.
* `costs`: the cost for each line
* `model`: a JuMP model
* `budget`: the budget for the budget constraint, which is fixed by `NDM.optimize()`
* `x`: x[l, f] indicates whether a line l is operated at frequency level f
* `θ`: θ[u, v] indicates the fraction of commuters for commute (u, v) can be serviced
* `solver`: GurobiSolver by default
* `options`: struct containing options for the optimization
"""
mutable struct MasterProblem
    # Network information
    np::TransitNetwork
    linelist::Vector{Vector{Int}}
    commutelines::Vector{Dict{Tuple{Int,Int},Any}}
    costs::Vector{Float64}
    
    # Optimization information
    model::JuMP.Model
    budget::JuMP.Variable
    x::Array{JuMP.Variable,2}
    θ::JuMP.JuMPArray
    solver

    # Parameters
    options::MasterOptions
end

"""
Returns a MasterProblem object for the given TransitNetwork, initialized on `linelist`
with default parameters in `options`.
To start with an empty set of lines, use linelist = Vector{Vector{Int}}()
"""
function MasterProblem(
        np::TransitNetwork;
        linelist::Vector{Vector{Int}} = uniquelines(np.lines),
        options = MasterOptions(),
        solver = Gurobi.GurobiSolver(OutputFlag = 0)
    )
    @assert in(options.modeltype, [:lp, :ip])
    @assert in(np.gridtype, [:latlong, :euclidean])
    @assert abs(options.angleparam) <= 1.0
    @assert options.distparam >= 1.0
    @assert length(options.freqwts) == length(options.xfrwts) == length(options.costwts) == options.nfreqs

    commutelines = allcommutelines(np, options.nlegs, linelist)
    costs = Vector{Float64}([linecost(np, line) for line in linelist])
    rmp, budget, x, θ = mastermodel(
        np, linelist, commutelines, costs, solver, options
    )
    
    MasterProblem(
        np, linelist, commutelines, costs, 
        rmp, budget, x, θ, solver, options
    )
end

"Builds base master problem model"
function mastermodel(
        np::TransitNetwork,
        linelist::Vector{Vector{Int}},
        commutelines::Vector{Dict{Tuple{Int,Int},Any}},
        costs::Vector{Float64},
        solver,
        options::MasterOptions
    )
    nlines = length(linelist)
    nfreqs = options.nfreqs
    freqwts = options.freqwts
    xfrwts = options.xfrwts
    costwts = options.costwts

    rmp = JuMP.Model(solver=solver)

    if options.modeltype == :lp
        JuMP.@variable(rmp, x[l=1:nlines,f=1:nfreqs] >= 0)
    elseif options.modeltype == :ip
        JuMP.@variable(rmp, x[l=1:nlines,f=1:nfreqs], Bin)
    end
    if options.nlegs == 2
        JuMP.@variable(rmp, twoline[
            (u,v)=commutes(np),
            w=np.xfrstns[(u,v)] 
        ])
    end
    JuMP.@variable(rmp, budget)
    JuMP.@variable(rmp, θ[(u,v)=commutes(np)])
    # maximize ridership
    JuMP.@objective(rmp, 
        Max, 
        sum(demand(np,(u,v))*θ[(u,v)] for (u,v) in commutes(np))
    )
    # choice constraints
    JuMP.@constraint(rmp,
        choseub[(u,v)=commutes(np)],
        θ[(u,v)] <= 1
    )
    if options.nlegs == 1
        JuMP.@constraint(rmp,
            choseline[(u,v)=commutes(np)],
            θ[(u,v)] <= sum(freqwts[f]*sum(x[l,f] for l in commutelines[1][u,v]) 
                            for f in 1:nfreqs))
    else 
        JuMP.@constraint(rmp,
            choseline[(u,v)=commutes(np)],
            θ[(u,v)] <= 
                sum(freqwts[f]*sum(x[l,f] for l in commutelines[1][u,v]) for f in 1:nfreqs) + 
                sum(twoline[(u,v),w] for w in np.xfrstns[(u,v)]))
        JuMP.@constraint(rmp,
            freq2a[(u,v)=commutes(np), w=np.xfrstns[(u,v)]],
            twoline[(u,v),w] <= 
                sum(xfrwts[f]*sum(x[l,f] for l in setdiff(commutelines[1][u,w], commutelines[1][u,v]))
                    for f in 1:nfreqs))
        JuMP.@constraint(rmp,
            freq2b[(u,v)=commutes(np), w=np.xfrstns[(u,v)]],
            twoline[(u,v),w] <= 
                sum(xfrwts[f]*sum(x[l,f] for l in setdiff(commutelines[1][w,v], commutelines[1][u,v]))
                    for f in 1:nfreqs))
    end

    if options.constrainedg
        edgelines = Dict{Tuple{Int,Int},Vector{Int}}()
        for l in 1:nlines, k in 2:length(linelist[l])
            u = linelist[l][k-1]
            v = linelist[l][k]
            edg = (min(u,v), max(u,v))
            if haskey(edgelines, edg)
                push!(edgelines[edg], l)
            else
                edgelines[edg] = [l]
            end
        end
        edges = keys(edgelines)
        coeffs = Dict{Tuple{Int,Int},Vector{Float64}}()
        for edg in edges
            if length(edgelines[edg]) == 1
                delete!(edgelines, edg)
            else
                coeffs[edg] = ones(length(edgelines[edg]))
                # allow for overlap between lines in original key network
                denom = length(intersect(1:options.ninitial, edgelines[edg]))
                if denom > 1
                    for l in 1:length(edgelines[edg])
                        if edgelines[edg][l] <= options.ninitial
                            coeffs[edg][l] = 1/denom
                        end
                    end
                end
            end
        end
        JuMP.@constraint(rmp,
            ccon[edg in keys(edgelines)], 
            sum(sum(coeffs[edg][l]*x[edgelines[edg][l],f] 
                    for l in 1:length(edgelines[edg])) 
                for f in 1:nfreqs) <= 1)
    end
    
    # choose one frequency
    if nfreqs > 1
        JuMP.@constraint(rmp, 
            [l=1:nlines],
            sum(x[l,f] for f in 1:nfreqs) <= 1)
    end

    # budget constraint
    JuMP.@constraint(rmp, bcon, 
        sum(costwts[f]*costs[l]*x[l,f] for l in 1:nlines for f in 1:nfreqs) <= budget)

    rmp, budget, x, θ
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
        rmp.np, rmp.linelist, rmp.commutelines, line
    )
    push!(rmp.costs, linecost(rmp.np, line))
    initialbudget = JuMP.getvalue(rmp.budget)

    rmp.model, rmp.budget, rmp.x, rmp.θ = 
        mastermodel(
            rmp.np, rmp.linelist, rmp.commutelines, rmp.costs, rmp.solver,
            rmp.options
        )
    JuMP.fix(rmp.budget, initialbudget)
end

"""
Returns a vector (with length `nlegs<=2`) of Dict{Tuple{Int,Int},Vector}, where

    commutelines[1][u,v] is a vector of the i lines that connect u to v, and
    commutelines[2][u,v] is empty (unused code)

"""
function allcommutelines(
        np::TransitNetwork,
        nlegs::Int,
        linelist::Vector{Vector{Int}}
    )
    @assert nlegs <= 2

    commutelines = [Dict{Tuple{Int,Int},Any}() for i in 1:nlegs]
    # 1. Initialize each (u,v) entry as an empty vector.
    for u=1:np.nstations, v=setdiff(1:np.nstations, u)
        commutelines[1][u,v] = Int[]
    end
    # 2. Populate the (u,v) entries in commutelines.
    # 
    # We begin with an empty list, and iteratively add each line using the
    # addline!() method. At the end, we should recover the original linelist.
    linelistcopy = Vector{Int}[]
    for line in linelist
        linelistcopy = addline!(
            np, linelistcopy, commutelines, line
        )
    end
    @assert linelistcopy == linelist

    return commutelines
end

"""
Modify `oldlines` and `commutelines` in-place, adding information about `line`.
"""
function addline!(
        np::TransitNetwork,
        oldlines::Vector{Vector{Int}},
        commutelines::Vector{Dict{Tuple{Int,Int},Any}},
        line::Vector{Int}
    )
    l1 = length(oldlines) + 1 # We introduce a new `line` with index `l1`.
    # Single-leg commutes
    for u in line, v in line
        u != v && push!(commutelines[1][u,v], l1)
    end
    # Two-leg commutes are now handed directly in the master problem, so no need to update.
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
