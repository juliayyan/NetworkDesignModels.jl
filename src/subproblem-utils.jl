"""
Computes all potential transfer stations between commutes (u,v).

It only deals with single-transfer commutes, and ignores all `(u,v)` pairs that
have a direct connection.

### Returns
A `(xfrstops_uw, xfrstops_wv)` tuple, where

* `xfrstops_uw[u,v]` is the vector of ints corresponding to `w`, and
* `xfrstops_wv[u,v]` is the vector of ints corresponding to `w`.
"""
function computexfrstns(rmp::MasterProblem, gridtype::Symbol)
    const linelist = rmp.linelist[find(round.(JuMP.getvalue(rmp.x), 5))]
    const nstns = rmp.np.nstations
    const stnlines = [find(in(u, line) for line in linelist) for u in 1:nstns]
    xfrstops_uw = Dict{Tuple{Int,Int},Vector{Int}}()
    xfrstops_wv = Dict{Tuple{Int,Int},Vector{Int}}()
    for u in 1:nstns, v in nonzerodests(rmp.np, u)
        xfrstops_uw[u,v] = Int[]
        xfrstops_wv[u,v] = Int[]
        # Ignore all (u,v) pairs that have a direct connection.
        length(intersect(stnlines[u], stnlines[v])) > 0 && continue
        # Deal with all (u,v) pairs with a single transfer station w.
        for w in 1:nstns
            if validtransfer(rmp.np, u, v, w, rmp.transferparam, gridtype)
                if length(intersect(stnlines[u], stnlines[w])) > 0
                    push!(xfrstops_uw[u,v], w)
                end
                if length(intersect(stnlines[w], stnlines[v])) > 0
                    push!(xfrstops_wv[u,v], w)
                end
            end 
        end
    end
    xfrstops_uw, xfrstops_wv
end

"finds path from source to sink using subproblem `sp`"
function getpath(sp::SubProblem)
    if round(JuMP.getobjectivevalue(sp.model)) > 0
        np = sp.np
        visited = falses(np.nstations)
        source_val = findfirst(round.(JuMP.getvalue(sp.src)))
        sink_val = findfirst(round.(JuMP.getvalue(sp.snk)))
        return getpath(
            source_val, source_val, sink_val, sp.edg, visited, sp.outneighbors
        )
    end
    Int[]
end

"helper function for wrapper getpath(), where `cur` can be 
 any starting node (in case of subtours)"
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
