"""
SubProblem
"""
mutable struct SubProblem
    np::TN.TransitNetworkProblem
    model::JuMP.Model
    src::Vector{JuMP.Variable}
    snk::Vector{JuMP.Variable}
    edg::JuMP.JuMPDict{JuMP.Variable}
    srv::JuMP.JuMPDict{JuMP.Variable}
    srv_uw
    srv_wv
    xfrstops_uw
    xfrstops_wv
    dists::Dict{Tuple{Int,Int},Float64}
    outneighbors::Vector{Vector{Int}}
    inneighbors::Vector{Vector{Int}}
    nlegs
    auxinfo::Dict{Symbol,Any}
end

"""
SubProblem Constructor
* Constructs an acyclic graph that only has edges 
* within some `delta` tolerance of `direction`
"""
function SubProblem(
    rmp::MasterProblem;
    nlegs::Int = length(rmp.commutelines), # want this in case 2-legs hard to solve
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    maxdist::Float64 = 0.5,
    direction::Vector{Float64} = [0.0,1.0],
    delta::Float64 = 1.0,
    maxlength::Int = 30
    )
    
    const np = rmp.np
    const nstns = np.nstations
    const gridtype = rmp.gridtype

    # construct graph and ensure there are no cycles
    dists = Dict{Tuple{Int,Int},Float64}() 
    outneighbors = [Int[] for u in 1:nstns]
    inneighbors  = [Int[] for u in 1:nstns]
    graph = LightGraphs.DiGraph(nstns)
    for u in 1:nstns, v in (u+1):nstns 
        d = edgecost(np, u, v, gridtype)
        b = dir(np,u,v,gridtype)
        sim = dot(direction,b)/norm(direction)/norm(b)
        if (d < maxdist) && (1-abs(sim) < delta)
            if sim >= 0
                dists[u,v] = d 
                push!(outneighbors[u], v)
                push!(inneighbors[v], u)
                LightGraphs.add_edge!(graph, (u,v))
            else
                dists[v,u] = d
                push!(outneighbors[v], u)
                push!(inneighbors[u], v)    
                LightGraphs.add_edge!(graph, (v,u))
            end
        end 
    end 
    @assert !LightGraphs.is_cyclic(graph)

    sp, src, snk, edg, srv, ingraph = basemodel(np, inneighbors, outneighbors, maxlength, solver)

    if nlegs == 2
        xfrstops_uw, xfrstops_wv = computexfrstns(np, rmp.linelist, rmp.transferparam, gridtype)
        srv_uw, srv_wv = transfermodel(np, sp, srv, ingraph, xfrstops_uw, xfrstops_wv)
    else 
        xfrstops_uw = xfrstops_wv = nothing
        srv_uw = srv_wv = nothing
    end

    SubProblem(np, 
        sp, src, snk, edg, srv, srv_uw, srv_wv,
        xfrstops_uw, xfrstops_wv,
        dists, outneighbors, inneighbors,
        nlegs, Dict{Symbol,Any}())
end 

"""
SubProblem Constructor
* Constructs a graph (not necessarily acyclic) 
* and solves by cutting planes
"""
function SubProblemCP(
    rmp::MasterProblem;
    nodeset::Vector{Int} = Vector(1:rmp.np.nstations),
    nlegs::Int = 1,
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    maxdist::Float64 = 0.5,
    maxlength::Int = 30
    )

    nlegs != 1 && warning("Only tested for nlegs = 1.")

    const np = rmp.np
    const nstns = np.nstations
    const gridtype = rmp.gridtype

    # construct graph and ensure there are no cycles
    dists = Dict{Tuple{Int,Int},Float64}() 
    outneighbors = [Int[] for u in 1:nstns]
    inneighbors  = [Int[] for u in 1:nstns]
    for u in 1:nstns
        for v in (u+1):nstns 
            !in(u, nodeset) && !in(v, nodeset) && continue
            d = edgecost(np, u, v, gridtype)
            if (d < maxdist)
                dists[u,v] = d 
                push!(outneighbors[u], v)
                push!(inneighbors[v], u)
                dists[v,u] = d
                push!(outneighbors[v], u)
                push!(inneighbors[u], v)    
            end 
        end
    end

    sp, src, snk, edg, srv, ingraph = basemodel(np, inneighbors, outneighbors, maxlength, solver)
    auxinfo = Dict{Symbol,Any}()
    auxinfo[:nlazy] = 0

    function removecycles(cb)
        visited = falses(np.nstations) # whether a node has been visited
        visited[setdiff(1:np.nstations,find(JuMP.getvalue(ingraph)))] = true
        source_val = sink_val = 0
        for u in 1:np.nstations 
            if round(JuMP.getvalue(src[u])) > 0.1
                source_val = u 
                break
            end
        end
        for u in 1:np.nstations 
            if round(JuMP.getvalue(snk[u])) > 0.1
                sink_val = u
                break
            end
        end
        @assert source_val > 0
        @assert sink_val > 0
        @assert source_val != sink_val
        simplepath = getpath(source_val, source_val, sink_val, edg, visited, outneighbors) # bus line
        while length(find(.!visited)) > 0 # search for subtours
            cyclenodes = getpath(findfirst(.!visited), source_val, sink_val, edg, visited, outneighbors)
            if length(cyclenodes) > 1
                expr = sum(sum(edg[u,v] for v in intersect(outneighbors[u], cyclenodes)) for u in cyclenodes)
                JuMP.@lazyconstraint(cb, expr <= length(cyclenodes) - 1)
                auxinfo[:nlazy] = auxinfo[:nlazy] + 1
            end
        end
    end
    JuMP.addlazycallback(sp, removecycles)
   
    SubProblem(np, 
        sp, src, snk, edg, srv, nothing, nothing,
        nothing, nothing,
        dists, outneighbors, inneighbors,
        nlegs, auxinfo) 
end
