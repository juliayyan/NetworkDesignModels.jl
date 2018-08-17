function getpath(sp::SubProblem)
    if round(JuMP.getobjectivevalue(sp.model)) > 0
        np = sp.np
        visited = falses(np.nstations)
        source_val = findfirst(round.(JuMP.getvalue(sp.src)))
        sink_val = findfirst(round.(JuMP.getvalue(sp.snk)))
        edge_val = JuMP.getvalue(sp.edg)
        return getpath(source_val, 
            source_val, sink_val,
            edge_val,
            visited, sp.outneighbors)
    end
    return Int[]
end

function getpath(cur::Int, 
    source_val::Int, sink_val::Int, 
    edge_val, 
    visited::BitArray,
    outneighbors::Vector{Vector{Int}};
    maxiterations = 1000)
    pathnodes = Int[]
    push!(pathnodes, cur)
    visited[cur] = true
    for i in 1:maxiterations
        for nxt in outneighbors[cur]
            if round(edge_val[cur,nxt]) > 0.1
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

