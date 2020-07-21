"""
Constructs a base direct-route subproblem.

Use transfermodel(...) to the additional variables and constraints needed for 
modelling single-transfers.
"""
function basemodel(
        rmp::MasterProblem,
        inneighbors::Vector{Vector{Int}},
        outneighbors::Vector{Vector{Int}},
        maxlength::Int, # maximum number of edges in a path
        linelist::Vector{Vector{Int}},
        solver
    )
    np = rmp.np
    nstns = np.nstations
    
    sp = JuMP.Model(solver=solver)

    JuMP.@variables sp begin
        src[u=1:nstns], Bin
        snk[u=1:nstns], Bin
        edg[u=1:nstns, v=outneighbors[u]], Bin
        0 <= srv[(u,v)=commutes(rmp)] <= 1
    end
    
    # path constraints
    JuMP.@constraints sp begin
        sum(edg) <= maxlength
        sum(src) == 1
        sum(snk) == 1
        # degree of nodes should be at most one
        # (handles an edge case that the subtour elimination doesn't
        #  catch, which is loops within a bus line)
        [u=1:nstns], src[u] + sum(edg[v,u] for v in inneighbors[u]) <= 1
        [u=1:nstns], snk[u] + sum(edg[u,v] for v in outneighbors[u]) <= 1
        [u=1:nstns], src[u] + sum(edg[v,u] for v in inneighbors[u]) ==
                     snk[u] + sum(edg[u,v] for v in outneighbors[u])
        # avoid singletons
        [u=1:nstns], src[u] + snk[u] <= 1
    end

    # demand service
    JuMP.@expression(sp,
        ingraph[u=1:nstns],
        src[u] + sum(edg[u2,u] for u2 in inneighbors[u]))
    JuMP.@constraints sp begin
        [(u,v)=commutes(rmp)], srv[(u,v)] <= ingraph[u]
        [(u,v)=commutes(rmp)], srv[(u,v)] <= ingraph[v]
    end

    # remove any lines that already exist
    function removeduplicates(cb)
        visited = falses(np.nstations) # whether a node has been visited
        visited[setdiff(1:np.nstations, findall(JuMP.getvalue(ingraph) .> 0))] .= true
        source_val = findfirst(round.(JuMP.getvalue.(src)) .> 0.1)
        sink_val = findfirst(round.(JuMP.getvalue.(snk)) .> 0.1)
        @assert source_val > 0
        @assert sink_val > 0
        @assert source_val != sink_val
        simplepath = getpath( # bus line
            source_val, source_val, sink_val, edg, visited, outneighbors
        )
        # cycles found
        if round(JuMP.getvalue(sum(edg)), digits=3) > length(simplepath)
            return
        end
        if linein(simplepath, linelist)
            for pathelim in [simplepath, reverse(simplepath)]
                try
                    inexpr = sum(edg[pathelim[k-1],pathelim[k]] for k in 2:length(pathelim))
                    npathedges = JuMP.getvalue(inexpr)
                    outexpr = (length(edg) - npathedges) - (sum(edg) - inexpr)
                    JuMP.@lazyconstraint(cb, inexpr + outexpr <= length(edg) - length(simplepath)/2)
                catch
                end
            end
        end
    end
    JuMP.addlazycallback(sp, removeduplicates)

    sp, src, snk, edg, srv, ingraph
end

"""
Adds transferring variables and constraints to subproblem model `sp`.

Use basemodel(...) for constructing a subproblem model.
"""
function transfermodel(
        rmp::MasterProblem,
        sp::JuMP.Model,
        srv,
        ingraph
    )
    np = rmp.np
    nstns = np.nstations

    JuMP.@variable(sp, 0 <= srv2a[(u,v)=commutes(rmp),w=np.xfrstns[(u,v)]] <= 1)
    JuMP.@variable(sp, 0 <= srv2b[(u,v)=commutes(rmp),w=np.xfrstns[(u,v)]] <= 1)

    JuMP.@constraints sp begin
        [(u,v)=commutes(rmp),w=np.xfrstns[(u,v)]], srv2a[(u,v),w] <= ingraph[u]
        [(u,v)=commutes(rmp),w=np.xfrstns[(u,v)]], srv2a[(u,v),w] <= ingraph[w]
        [(u,v)=commutes(rmp),w=np.xfrstns[(u,v)]], srv2a[(u,v),w] <= 1-ingraph[v]
        [(u,v)=commutes(rmp),w=np.xfrstns[(u,v)]], srv2b[(u,v),w] <= ingraph[v]
        [(u,v)=commutes(rmp),w=np.xfrstns[(u,v)]], srv2b[(u,v),w] <= ingraph[w]
        [(u,v)=commutes(rmp),w=np.xfrstns[(u,v)]], srv2b[(u,v),w] <= 1-ingraph[u]
        [(u,v)=commutes(rmp)], srv[(u,v)] + sum(srv2a[(u,v),w] + srv2b[(u,v),w] for w in np.xfrstns[(u,v)]) <= 1
    end

    (srv2a, srv2b)
end

"""
Uses dual values `p` and `q` to generate a profitable line.

This function helps solve the subproblem many times just changing the objective
without rebuilding the model.

### Keyword Arguments
* `sp`: the SubProblem instance
* `rmp`: the MasterProblem from which we get the duals
* `f`: which frequency level to optimize for
* `trackingstatuses`: for tracking information through solver callbacks.
* `trackingtimegrid`: how often to save tracking data in seconds

### Returns
A `path::Vector{Int}` of the stations along the profitable line.
"""
function generatecolumn(
        sp::SubProblem, 
        rmp::MasterProblem;
        f::Int = 1,
        trackingstatuses::Vector{Symbol} = Symbol[],
        trackingtimegrid::Int = 5
    )
    # set objective of subproblem based on dual values
    freqwt = rmp.options.freqwts[f]
    xfrwt = rmp.options.xfrwts[f]
    costwt = rmp.options.costwts[f]
    p = JuMP.getdual(rmp.model[:choseline])
    q = max(1e-3,JuMP.getdual(rmp.model[:bcon]))
    edgexpr = 0
    if rmp.options.constrainedg
        for key in keys(rmp.model[:ccon])
            (u,v) = key[1]
            if in(v, sp.outneighbors[u])
                edgexpr += JuMP.getdual(rmp.model[:ccon][(u,v)]) * sp.edg[u,v]
            end
            if in(u, sp.outneighbors[v])
                edgexpr += JuMP.getdual(rmp.model[:ccon][(u,v)]) * sp.edg[v,u]
            end
        end
    end
    if sp.srv2 != nothing
        pi2a = JuMP.getdual(rmp.model[:freq2a])
        pi2b = JuMP.getdual(rmp.model[:freq2b])
        JuMP.@objective(sp.model,
            Max,
            freqwt * sum(p[(u,v)] * sp.srv[(u,v)] for (u,v) in commutes(rmp)) +
            xfrwt * sum(sum(pi2a[(u,v),w] * sp.srv2[1][(u,v),w] 
                    for w in rmp.np.xfrstns[(u,v)]) 
                for (u,v) in commutes(rmp)) +
            xfrwt * sum(sum(pi2b[(u,v),w] * sp.srv2[2][(u,v),w] 
                    for w in rmp.np.xfrstns[(u,v)]) 
                for (u,v) in commutes(rmp)) - 
            costwt * q * sum(sp.dists[u,v]*sp.edg[u,v]
                for u in 1:sp.np.nstations, v in sp.outneighbors[u]) - edgexpr
        )
    else
        JuMP.@objective(sp.model,
            Max,
            freqwt * sum(p[(u,v)] * sp.srv[(u,v)] for (u,v) in commutes(rmp)) -
            costwt * q * sum(sp.dists[u,v]*sp.edg[u,v]
                for u in 1:sp.np.nstations, v in sp.outneighbors[u]) - edgexpr
        )
    end

    # track solve information
    t0 = time()
    for tracking in trackingstatuses
        if tracking == :Intermediate
            sp.auxinfo[:obj]       = Float64[]
            sp.auxinfo[:bestbound] = Float64[]
            sp.auxinfo[:time]      = Float64[]
            function boundscallback(cb)
                currtime = time() - t0
                if ((length(sp.auxinfo[:time]) == 0) ||
                    (currtime > sp.auxinfo[:time][end] + trackingtimegrid))
                    push!(sp.auxinfo[:obj], MathProgBase.cbgetobj(cb))
                    push!(sp.auxinfo[:bestbound],
                          MathProgBase.cbgetbestbound(cb))
                    push!(sp.auxinfo[:time], time() - t0)
                end
            end
            JuMP.addinfocallback(sp.model, boundscallback, when = :Intermediate)
        elseif tracking == :MIPSol
            sp.auxinfo[:solntime] = Float64[]
            function solncallback(cb)
                push!(sp.auxinfo[:solntime], time() - t0)
            end
            JuMP.addinfocallback(sp.model, solncallback, when = :MIPSol)
        end
    end

    # solve
    JuMP.solve(sp.model)
    sp.auxinfo[:endtime] = time() - t0

    # save data
    if JuMP.getobjectivevalue(sp.model) > 0
        path = getpath(sp) 
        if ((length(path) > 0) &&
            (round(sum(JuMP.getvalue(sp.edg))) != length(path) - 1))
            error("Dual solution is not a valid path")
        end
    else
        path = Vector{Int}()
    end

    path
end 

"""
Heuristically generates columns.

It starts with edge-restricted networks and iteratively builds up
relevant sections of the network.

### Keyword Arguments
* `rmp`: the MasterProblem to generate a column for.  Builds subproblems internally
* `directions`: A vector of directions z to run the subproblem with
    preprocessed edge set E(z, delta).
* `nlegs`: The (maximum) number of legs of each commute.
* `maxdist`: The threshold distance for two stations to be considered
    neighbors. The distance is measured based on the `gridtype` in
    `rmp::MasterProblem`.
* `delta`: Tolerance within `direction` for edges in the graph to obey.
* `maxlength`: The maximum number of edges in a path.
* `stepsize`: The number of neighbor-lookaheads. `stepsize = 1` means that at
    each iteration, we build up the graph by looking at the current nodeset +
    their neighbors. `stepsize = 2` means that we look at the current nodeset +
    neighbors + neighbors' neighbors, etc.
* `solver`: The solver being used to solve the problem.
* `maxiterations`: The maximum number of iterations to build up relevant
        sections of the network for.
* `trackingstatuses`: For tracking information through solver callbacks.

### Returns
A `(path, auxinfo)` tuple.
"""
function generatecolumn(
        rmp::MasterProblem; 
        directions::Vector{Vector{Float64}} = [
            [ 1.0, 0.0], [ 1.0, 0.5],
            [ 1.0, 1.0], [ 0.5, 1.0],
            [ 0.0, 1.0], [-0.5, 1.0],
            [-1.0, 1.0], [-1.0,-0.5]
        ],
        nlegs::Int = length(rmp.commutelines),
        maxdist::Float64 = 0.5,
        delta::Float64 = 1.0,
        maxlength::Int = 30,
        solver = Gurobi.GurobiSolver(OutputFlag = 0),
        stepsize::Int = 1,
        maxiterations::Int = 10,
        trackingstatuses::Vector{Symbol} = Symbol[]
    )
    auxinfo = Dict{Symbol,Any}()
    auxinfo[:time]   = Float64[]
    auxinfo[:obj]    = Float64[]
    auxinfo[:nnodes] = Int[]
    auxinfo[:nlazy]  = 0
    t0 = time()

    nstns = rmp.np.nstations
    dists = [
        edgecost(rmp.np,u,v) for u in 1:nstns, v in 1:nstns
    ]
    neighbors = [setdiff(find(dists[u,:] .< maxdist),u) for u in 1:nstns]

    # compute warm start
    oldobj = 0.0
    sp_warm = [SubProblem(rmp,
                          nlegs = nlegs,
                          maxdist = maxdist, 
                          maxlength = maxlength,
                          direction = d, 
                          solver = Gurobi.GurobiSolver(OutputFlag = 0)) 
                for d in directions]
    sp_objs = [JuMP.getobjectivevalue(spw.model) for spw in sp_warm]
    soln_warm = [generatecolumn(spw,rmp) for spw in sp_warm]
    nodeset = unique(union(soln_warm...))
    push!(auxinfo[:time], time() - t0)
    push!(auxinfo[:obj], maximum(sp_objs))
    push!(auxinfo[:nnodes], length(nodeset))

    # iteratively generate path
    path = soln_warm[findmax(sp_objs)[2]]
    for i in 1:maxiterations
        sp = SubProblemCP(rmp, 
                          nlegs = nlegs,
                          maxdist = maxdist, 
                          maxlength = maxlength,
                          nodeset = nodeset, 
                          solver = solver)
        warmstart(sp, path)
        path = generatecolumn(sp,rmp,
                              trackingstatuses=trackingstatuses)
        if round(JuMP.getobjectivevalue(sp.model),digits=3) == round(oldobj,digits=3)
            break
        else
            auxinfo[:nlazy] += sp.auxinfo[:nlazy]
            push!(auxinfo[:time], time() - t0)
            push!(auxinfo[:obj],  JuMP.getobjectivevalue(sp.model))
            nodeset = union(nodeset, path)
            for j in 2:stepsize 
                nodeset = union(nodeset, neighbors[nodeset]...)
            end
            push!(auxinfo[:nnodes], length(nodeset))
            oldobj = JuMP.getobjectivevalue(sp.model)
        end
    end

    auxinfo[:endtime] = time() - t0
    
    path, auxinfo 
end

"Warm starts a SubProblem with a given path"
function warmstart(sp::SubProblem, path::Vector{Int})
    JuMP.setvalue(sp.src[path[1]]  , 1)
    JuMP.setvalue(sp.snk[path[end]], 1)
    for i in 2:length(path)
        JuMP.setvalue(sp.edg[path[i-1],path[i]], 1)
    end
end

"""
Warm starts a SubProblem with a proposed path from which we can construct a partial solution
Necessary in case the proposed path is not actually feasible for the SubProblem
"""
function trywarmstart(sp::SubProblem, trypath::Vector{Int})
    warmstarts = Vector{Vector{Int}}()
    curpath = Vector{Int}()
    for k in 2:length(trypath)
        if in(trypath[k], sp.outneighbors[trypath[k-1]])
            if length(curpath) == 0
                push!(curpath, trypath[k-1])
            end
            push!(curpath, trypath[k])
        elseif length(curpath) > 0
            push!(warmstarts, curpath)
            curpath = Vector{Int}()
        end
        if (k == length(trypath)) && (length(curpath) > 0)
            push!(warmstarts, curpath)
        end
    end
    objs = Vector{Float64}()
    for path in warmstarts
        obj = 0
        for u in path, v in path
            try
                obj += JuMP.getdual(rmp.model[:choseline][(u,v)])
            catch
            end
        end
        push!(objs, obj)
    end
    if length(objs) > 0
        warmstart(sp, warmstarts[findmax(objs)[2]])
    end
end
