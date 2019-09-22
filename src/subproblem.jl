mutable struct SubProblem
    np::TN.TransitNetworkProblem
    model::JuMP.Model
    src::Vector{JuMP.Variable}
    snk::Vector{JuMP.Variable}
    edg::JuMP.JuMPDict{JuMP.Variable}
    ingraph
    srv::JuMP.JuMPDict{JuMP.Variable}
    srv2
    xfrstops_uw
    xfrstops_wv
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
* `nlegs`: the (maximum) number of legs of each commute.
* `solver`: The solver being used to solve the problem.
* `maxdist`: the threshold distance for two stations to be considered.
    neighbors. The distance is measured based on the `gridtype` in
    `rmp::MasterProblem`.
* `direction`: The direction for which we run the subproblem with preprocessed
    edge set E(direction, delta).
* `delta`: tolerance within `direction` for edges in the graph to obey.
* `maxlength`: the maximum number of edges in a path.
"""
function SubProblem(
        rmp::MasterProblem;
        nlegs::Int = length(rmp.commutelines), # in case 2-legs hard to solve
        solver = Gurobi.GurobiSolver(OutputFlag = 0),
        maxdist::Float64 = 0.5,
        direction::Vector{Float64} = [0.0,1.0],
        delta::Float64 = 1.0,
        maxlength::Int = 30 # maximum number of edges in a path
    )
    np = rmp.np
    nstns = np.nstations
    gridtype = rmp.gridtype

    # construct graph and ensure there are no cycles
    dists = Dict{Tuple{Int,Int},Float64}()
    outneighbors = [Int[] for u in 1:nstns]
    inneighbors  = [Int[] for u in 1:nstns]
    graph = LightGraphs.DiGraph(nstns)
    for u in 1:nstns, v in (u+1):nstns
        d = edgecost(np, u, v, gridtype)
        b = dir(np, u, v, gridtype)
        sim = dot(direction, b) / norm(direction) / norm(b)
        if (d < maxdist) && (1 - abs(sim) < delta)
            if sim >= 0
                dists[u,v] = d
                push!(outneighbors[u], v)
                push!(inneighbors[v], u)
                LightGraphs.add_edge!(graph, (u, v))
            else
                dists[v,u] = d
                push!(outneighbors[v], u)
                push!(inneighbors[u], v)
                LightGraphs.add_edge!(graph, (v, u))
            end
        end
    end
    @assert !LightGraphs.is_cyclic(graph)

    sp, src, snk, edg, srv, ingraph = basemodel(
        np, inneighbors, outneighbors, maxlength, solver
    )

    if nlegs == 2
        xfrstops_uw, xfrstops_wv = computexfrstns(rmp, gridtype)
        srv2 = transfermodel(
            np, sp, srv, ingraph, xfrstops_uw, xfrstops_wv
        )
    else
        srv2 = xfrstops_uw = xfrstops_wv = nothing
    end

    SubProblem(
        np, sp, src, snk, edg, ingraph, srv, srv2, xfrstops_uw, xfrstops_wv,
        dists, outneighbors, inneighbors, nlegs, Dict{Symbol,Any}()
    )
end

"""
SubProblem with Cutting Planes

Constructs a graph (not necessarily acyclic) to be optimized over by using
cutting planes.

### Keyword Arguments
* `nodeset`: The set of stations under consideration.
* `nlegs`: the (maximum) number of legs of each commute.
* `solver`: The solver being used to solve the problem.
* `maxdist`: the threshold distance for two stations to be considered.
    neighbors. The distance is measured based on the `gridtype` in
    `rmp::MasterProblem`.
* `maxlength`: the maximum number of edges in a path.
"""
function SubProblemCP(
        rmp::MasterProblem;
        nodeset::Vector{Int} = Vector(1:rmp.np.nstations),
        nlegs::Int = length(rmp.commutelines),
        solver = Gurobi.GurobiSolver(OutputFlag = 0),
        maxdist::Float64 = 0.5,
        maxlength::Int = 30,
        traveltimes::Bool = false,
        detail::Bool = false
    )
    np = rmp.np
    nstns = np.nstations
    gridtype = rmp.gridtype

    # construct graph and ensure there are no cycles
    dists = Dict{Tuple{Int,Int},Float64}() 
    outneighbors = [Int[] for u in 1:nstns]
    inneighbors  = [Int[] for u in 1:nstns]
    for u in 1:nstns
        for v in (u+1):nstns 
            !in(u, nodeset) && !in(v, nodeset) && continue
            d = edgecost(np, u, v, gridtype)
            if d < maxdist
                dists[u,v] = d
                push!(outneighbors[u], v)
                push!(inneighbors[v], u)
                dists[v,u] = d
                push!(outneighbors[v], u)
                push!(inneighbors[u], v)
            end
        end
    end

    sp, src, snk, edg, srv, ingraph = basemodel(
        np, inneighbors, outneighbors, maxlength, solver
    )
    if nlegs == 2
        xfrstops_uw, xfrstops_wv = computexfrstns(rmp, gridtype)
        srv2 = transfermodel(
            np, sp, srv, ingraph, xfrstops_uw, xfrstops_wv
        )
    else
        srv2 = xfrstops_uw = xfrstops_wv = nothing
    end

    auxinfo = Dict{Symbol,Any}()
    auxinfo[:nlazy] = 0

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
                expr = sum(sum(edg[u,v]
                               for v in intersect(outneighbors[u], cyclenodes))
                           for u in cyclenodes)
                JuMP.@lazyconstraint(cb, expr <= length(cyclenodes) - 1)
                auxinfo[:nlazy] += 1
            end
        end
        if traveltimes
            addlazytraveltimes(cb, 
                np, edg, outneighbors,
                rmp.gridtype, rmp.distparam,
                detail,
                auxinfo,
                simplepath)
        end
    end
    JuMP.addlazycallback(sp, removecycles)
   
    SubProblem(
        np, sp, src, snk, edg, ingraph, srv, srv2, xfrstops_uw, xfrstops_wv,
        dists, outneighbors, inneighbors, nlegs, auxinfo
    )
end
