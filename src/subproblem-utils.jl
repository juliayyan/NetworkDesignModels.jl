"""
Returns path from source to sink.

It assumes that the subproblem (being passed in) has been solved to a solution
with a positive objective value. Returns an empty vector otherwise.

### Arguments
* `sp`: SubProblem whose optimal solution we're retrieving as a path.

### Returns
A vector of ints corresponding to stations along the path.
"""
function getpath(sp::SubProblem)
    if round(JuMP.getobjectivevalue(sp.model)) > 0
        visited = falses(sp.np.nstations)
        source_val = findfirst(round.(JuMP.getvalue(sp.src)) .> 0.1)
        sink_val = findfirst(round.(JuMP.getvalue(sp.snk)) .> 0.1)
        return getpath(
            source_val, source_val, sink_val, sp.edg, visited, sp.outneighbors
        )
    end

    Int[]
end

"""
Helper function for getpath(), where `cur` can be any starting node (in case of
subtours).
"""
function getpath(
        cur::Int, 
        source_val::Int,
        sink_val::Int, 
        edg::JuMP.JuMPDict{JuMP.Variable}, 
        visited::BitArray,
        outneighbors::Vector{Vector{Int}};
        maxiterations = 1000
    )
    pathnodes = Int[]
    push!(pathnodes, cur)
    visited[cur] = true
    for i in 1:maxiterations
        for nxt in outneighbors[cur]
            if round(JuMP.getvalue(edg[cur,nxt])) > 0.1
                if in(nxt, pathnodes)
                    return pathnodes
                elseif nxt == sink_val
                    push!(pathnodes, nxt)
                    visited[nxt] = true
                    return pathnodes
                end
                push!(pathnodes, nxt)
                visited[nxt] = true
                cur = nxt
                break
            end 
        end
    end
    
    pathnodes
end
