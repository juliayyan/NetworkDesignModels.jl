"""
Calculates traveltimes for commutes on the line and generates cuts if any commutes
have travel times that are too long.
"""
function addlazytraveltimes(
    cb,
    np::TransitNetwork,
    ingraph,
    edg,
    outneighbors::Vector{Vector{Int}},
    distparam::Float64,
    detail::Bool,
    auxinfo::Dict{Symbol,Any},
    simplepath::Vector{Int})
    for len in 2:length(simplepath), p_i in 1:(length(simplepath)-len)
        p_j = p_i+len
        directdist = spcost(np,simplepath[p_i],simplepath[p_j]) 
        pathsubset = simplepath[p_i:p_j]
        subsetcost = linecost(np,pathsubset)
        if subsetcost > distparam*directdist
            # calculate the shortest path connecting all nodes in this subset
            if detail
                mincost = minimum(insertionheuristic(np,pathsubset,start)[1]
                                  for start in 1:length(pathsubset))
            else 
                mincost = insertionheuristic(np,pathsubset)[1]
            end
            if mincost >= distparam*directdist
                # cut all paths containing all of these nodes
                expr = sum(ingraph[u] for u in pathsubset)
                auxinfo[:travelshort] += 1
            else
                expr = 0
                for j in 1:length(pathsubset)-1, k in (j+1):length(pathsubset)
                    if in(pathsubset[k], outneighbors[pathsubset[j]])
                        expr += edg[pathsubset[j], pathsubset[k]]
                    end
                end
                auxinfo[:traveltourn] += 1
            end
            JuMP.@lazyconstraint(cb, 
                expr <= JuMP.getvalue(expr) - 1)
            auxinfo[:nlazy] += 1
            break
        end
    end
end

"""
uses an insertion heuristic to calculate shortest path through nodes
starting at nodes[start]
"""
function insertionheuristic(
    np::TransitNetwork, 
    nodes::Vector{Int},
    start::Int = 1)
    nnodes = length(nodes)
    @assert length(unique(nodes)) == nnodes
    visited = falses(nnodes)
    totalcost = 0
    this = start
    visited[start] = true
    path = Vector{Int}()
    push!(path, start)
    while sum(.!visited) > 0
        costs = [visited[that] ? Inf : spcost(np,nodes[this],nodes[that]) 
                 for that in 1:nnodes]
        mincost, that = findmin(costs)
        if mincost == Inf
            return Inf, Int[]
        end
        totalcost += mincost
        this = that
        push!(path, this)
        visited[that] = true
    end
    totalcost, nodes[path]
end

"""
Returns arrays of commutes that cannot all coexist due to travel time restrictions
"""
function badcommutecombos(np::TransitNetwork, options::MasterOptions; 
    rmp = nothing,
    k::Int = 2,
    lb::Int = 100,
    prevnodes::Vector{Vector{Int}} =  Vector{Vector{Int}}())
    @assert k <= 4
    coms = commutes(np, lb)
    if rmp != nothing
        coms = intersect(coms, commutes(rmp))
    end
    badcoms = []
    indices = IterTools.subsets(1:length(coms), k)
    ProgressMeter.@showprogress for ids in indices
        nodes = unique(vcat([vcat(com...) for com in coms[ids]]...))
        if any([issubset(nodes, pn) for pn in prevnodes])
            continue
        end
        !checknodes(np, options.distparam, nodes) && push!(badcoms, coms[ids])
    end
    badcoms
end

"check whether some path connecting all of `nodes` is feasible for travel time restrictions"
function checknodes(np::TransitNetwork, distparam::Float64, nodes::Vector{Int})
    nnodes = length(nodes)
    if nnodes <= 4
        nodepaths = Combinatorics.permutations(nodes)
    else
        tries = [insertionheuristic(np,nodes,k) for k in 1:nnodes]
        costs = [tr[1] for tr in tries]
        nodepaths = [tr[2] for tr in tries]
        nodepaths = [nodepaths[findmin(costs)[2]]]
    end
    for path in nodepaths
        if checkline(np, distparam, path)
            return true
        end
    end
    false
end

"check whether the given `path` is feasible for travel time restrictions"
function checkline(np::TransitNetwork, distparam::Float64, path::Vector{Int})
    nnodes = length(path)
    for o in 1:nnodes, d in (o+1):nnodes
        thiscost = sum([spcost(np, path[l-1], path[l]) for l in (o+1):d])
        mincost = spcost(np, path[o], path[d])
        if thiscost > distparam * mincost
            return false
        end
    end
    true
end
