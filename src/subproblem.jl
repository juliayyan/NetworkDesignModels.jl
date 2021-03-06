"""
SubProblem

* `np`: a TransitNetwork object that contains the network information
* `model`: a JuMP model
* `src`: src[u] indicates whether the source is at stop u
* `snk`: snk[u] indicates whether the sink is at stop v
* `edg`: edg[u,v] indicates whether edge (u,v) is used
* `ingraph`: ingraph[u] is an expression indicating whether u is used in the line
* `srv`: srv[u,v] indicates whether commute (u,v) can be served by the line directly (one leg)
* `srv2`: srv2[u,v] indicates whether commute (u,v) can be served by the line with one transfer (two legs)
* `dists`: dists[u,v] is the distance of edge (u,v)
* `outneighbors`: outneighbors[u] is a vector containing other stations v for which (u, v) is a valid edge
* `inneighbors`: inneighbors[u] is a vector containing other stations v for which (v, u) is a valid edge
* `nlegs`: passed by the MasterProblem; 1 (direct-route model) or 2 (single-transfer model)
* `auxinfo`: stores auxiliary solver information
"""
mutable struct SubProblem
    np::TransitNetwork
    model::JuMP.Model
    src::Vector{JuMP.Variable}
    snk::Vector{JuMP.Variable}
    edg::JuMP.JuMPDict{JuMP.Variable}
    ingraph
    srv::JuMP.JuMPArray
    srv2
    dists::Dict{Tuple{Int,Int},Float64}
    outneighbors::Vector{Vector{Int}}
    inneighbors::Vector{Vector{Int}}
    nlegs::Int
    auxinfo::Dict{Symbol,Any}
end

"""
SubProblem with Edge Processing

Constructs an acyclic graph that only has edges within some `delta` tolerance of
`direction`.

### Keyword Arguments
* `nlegs`: the (maximum) number of legs of each commute (should be same as in the MasterProblem).
* `solver`: The solver being used to solve the problem.
* `maxdist`: the threshold distance for two stations to be considered.
    neighbors. The distance is measured based on the `gridtype` in
    `rmp::MasterProblem`.
* `direction`: The direction for which we run the subproblem with preprocessed
    edge set E(direction, delta).
* `delta`: tolerance within `direction` for edges in the graph to obey.
* `maxlength`: the maximum number of edges in a path.
* `traveltimes`: whether to constrain travel times on the network
* `detail`: whether to use a more detailed heuristic to compute shortest paths for travel times
"""
function SubProblem(
        rmp::MasterProblem;
        nlegs::Int = rmp.options.nlegs,
        solver = Gurobi.GurobiSolver(OutputFlag = 0),
        maxdist::Float64 = 0.5,
        direction::Vector{Float64} = [0.0,1.0],
        delta::Float64 = 1.0, # 1 - (delta in paper) -- 1.0 in the code corresponds to 0.0 in the paper
        maxlength::Int = 30, # maximum number of edges in a path
        traveltimes::Bool = false,
        detail::Bool = false
    )
    np = rmp.np
    nstns = np.nstations

    # construct graph and ensure there are no cycles
    if rmp.np.dists == nothing
        dists = Dict{Tuple{Int,Int},Float64}()
    else
        dists = rmp.np.dists
    end
    outneighbors = [Int[] for u in 1:nstns]
    inneighbors  = [Int[] for u in 1:nstns]
    graph = LightGraphs.DiGraph(nstns)
    for u in 1:nstns, v in (u+1):nstns
        d = edgecost(np, u, v)
        b = dir(np, u, v)
        sim = dot(direction, b) / norm(direction) / norm(b)
        if (d < maxdist) && (1 - abs(sim) <= delta)
            if sim >= 0
                if rmp.np.dists == nothing
                    dists[u,v] = d
                end
                push!(outneighbors[u], v)
                push!(inneighbors[v], u)
                LightGraphs.add_edge!(graph, (u, v))
            else
                if rmp.np.dists == nothing
                    dists[v,u] = d
                end
                push!(outneighbors[v], u)
                push!(inneighbors[u], v)
                LightGraphs.add_edge!(graph, (v, u))
            end
        end
    end
    @assert !LightGraphs.is_cyclic(graph)

    # construct base JuMP model
    sp, src, snk, edg, srv, ingraph = basemodel(
        rmp, inneighbors, outneighbors, maxlength, rmp.linelist, solver
    )
    # add variables and constraints for single-transfer model if relevant
    if nlegs == 2
        srv2 = transfermodel(rmp, sp, srv, ingraph)
    else
        srv2 = nothing
    end

    # initializes empty dictionary for auxiliary solve info to be stored
    auxinfo = initializedict()

    # adds cuts for limiting travel time on network
    function cuttraveltimes(cb)
        # get line from the variable values
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
        # add the cut
        if traveltimes
            addlazytraveltimes(cb, 
                np, ingraph, edg, outneighbors, rmp.options.distparam,
                detail,
                auxinfo,
                simplepath)
        end
    end
    if traveltimes
        JuMP.addlazycallback(sp, cuttraveltimes)
    end

    SubProblem(
        np, sp, src, snk, edg, ingraph, srv, srv2,
        dists, outneighbors, inneighbors, nlegs, auxinfo
    )
end

"""
SubProblem with Cutting Planes

Constructs a graph (not necessarily acyclic) to be optimized over by using
cutting planes to remove subtours.

### Keyword Arguments
* `nodeset`: The set of stations under consideration (default all).
* `nlegs`: the (maximum) number of legs of each commute (should be same as in the MasterProblem).
* `solver`: The solver being used to solve the problem.
* `maxdist`: the threshold distance for two stations to be considered.
    neighbors. The distance is measured based on the `gridtype` in
    `rmp::MasterProblem`.
* `maxlength`: the maximum number of edges in a path.
* `traveltimes`: whether to constrain travel times on the network
* `detail`: whether to use a more detailed heuristic to compute shortest paths for travel times
"""
function SubProblemCP(
        rmp::MasterProblem;
        nodeset::Vector{Int} = Vector(1:rmp.np.nstations),
        nlegs::Int = rmp.options.nlegs,
        solver = Gurobi.GurobiSolver(OutputFlag = 0),
        maxdist::Float64 = 0.5,
        maxlength::Int = 30,
        traveltimes::Bool = false,
        detail::Bool = false
    )
    np = rmp.np
    nstns = np.nstations

    # construct graph and ensure there are no cycles
    if rmp.np.dists == nothing
        dists = Dict{Tuple{Int,Int},Float64}()
    else
        dists = rmp.np.dists
    end
    outneighbors = [Int[] for u in 1:nstns]
    inneighbors  = [Int[] for u in 1:nstns]
    for u in 1:nstns
        for v in (u+1):nstns 
            !in(u, nodeset) && !in(v, nodeset) && continue
            d = edgecost(np, u, v)
            if d < maxdist
                if rmp.np.dists == nothing
                    dists[u,v] = d
                    dists[v,u] = d
                end
                push!(outneighbors[u], v)
                push!(inneighbors[v], u)
                push!(outneighbors[v], u)
                push!(inneighbors[u], v)
            end
        end
    end

    # construct base JuMP model
    sp, src, snk, edg, srv, ingraph = basemodel(
        rmp, inneighbors, outneighbors, maxlength, rmp.linelist, solver
    )
    # add variables and constraints for single-transfer model if relevant
    if nlegs == 2
        srv2 = transfermodel(rmp, sp, srv, ingraph)
    else
        srv2 = nothing
    end
    # some easy cuts to start
    JuMP.@constraint(sp, 
        [u=1:np.nstations,v=outneighbors[u]],
        edg[u,v] + edg[v,u] <= 1)

    # initializes empty dictionary for auxiliary solve info to be stored
    auxinfo = initializedict()
    
    # subtour elimination
    function removecycles(cb)
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
        while length(findall(.!visited)) > 0 # search for subtours
            cyclenodes = getpath(
                findfirst(.!visited), source_val, sink_val, edg, visited,
                outneighbors
            )
            if length(cyclenodes) > 1
                if traveltimes && !checknodes(np, rmp.options.distparam, cyclenodes)
                    expr = sum(ingraph[u] for u in cyclenodes)
                else
                    expr = sum(sum(edg[u,v]
                                   for v in intersect(outneighbors[u], cyclenodes))
                               for u in cyclenodes)
                end
                JuMP.@lazyconstraint(cb, expr <= length(cyclenodes) - 1)
                auxinfo[:cycle] += 1
                auxinfo[:nlazy] += 1
                break
            end
        end
        if traveltimes
            addlazytraveltimes(cb, 
                np, ingraph, edg, outneighbors, rmp.options.distparam,
                detail,
                auxinfo,
                simplepath)
        end
    end
    JuMP.addlazycallback(sp, removecycles)
   
    SubProblem(
        np, sp, src, snk, edg, ingraph, srv, srv2,
        dists, outneighbors, inneighbors, nlegs, auxinfo
    )
end

"""
Creates an empty Dict saving solve information
"""
function initializedict()
    auxinfo = Dict{Symbol,Any}()
    auxinfo[:nlazy] = 0
    auxinfo[:travelshort] = 0
    auxinfo[:traveltourn] = 0
    auxinfo[:cycle] = 0
    auxinfo[:checknodes] = Dict{Vector{Int},Bool}()
    auxinfo
end
