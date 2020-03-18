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
    end
    totalcost, path
end
