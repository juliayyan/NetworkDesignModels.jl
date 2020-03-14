function addlazytraveltimes(
    cb,
    np::TransitNetwork,
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
            if detail && length(pathsubset) <= 4
                mincost = minimum(insertionheuristic(np,pathsubset,start) 
                                  for start in 1:length(pathsubset))
            else 
                mincost = insertionheuristic(np,pathsubset)    
            end
            if mincost >= distparam*directdist
                expr = sum(
                    length(intersect(outneighbors[u], pathsubset)) == 0 ? 
                        0 : 
                        sum(edg[u,v] 
                            for v in intersect(outneighbors[u], pathsubset)) 
                        for u in pathsubset)
            else
                expr = 0
                for j in 1:length(pathsubset)-1, k in (j+1):length(pathsubset)
                    if in(pathsubset[k], outneighbors[pathsubset[j]])
                        expr += edg[pathsubset[j], pathsubset[k]]
                    end
                end
            end
            JuMP.@lazyconstraint(cb, 
                expr <= length(pathsubset) - 2)
            # this is wrong, too restrictive
            # sum(ingraph[u] for u in pathsubset) <= length(pathsubset) - 1)
            auxinfo[:nlazy] += 1
            break
        end
    end
end

"""
uses an insertion heuristic to calculate shortest path through nodes
starting at nodes[1]
"""
function insertionheuristic(
    np::TransitNetwork, 
    nodes::Vector{Int},
    start::Int = 1)
    @assert length(unique(nodes)) == length(nodes)
    visited = [false for u in nodes]
    totalcost = 0
    this = nodes[start]
    visited[start] = true
    while sum(.!visited) > 0
        candidates = nodes[.!visited]
        costs = [NetworkDesignModels.edgecost(np,this,that) 
                 for that in candidates]
        totalcost += minimum(costs)
        that = candidates[findmin(costs)[2]]
        this = that
        visited[findfirst(nodes .== that)] .= true
    end
    totalcost
end
