"""
Computes all potential transfer stations between commutes (u,v).

It only deals with single-transfer commutes, and ignores all `(u,v)` pairs that
have a direct connection.

### Returns
A `(xfrstops_uw, xfrstops_wv)` tuple, where

* `xfrstops_uw[u,v]` contains the `w`s that have a connection from `u`, 
    separated into vectors corresponding to active [1] and inactive [2] lines
* `xfrstops_wv[u,v]` contains the `w`s that have a connection to `v`,
    separated into vectors corresponding to active [1] and inactive [2] lines
"""
function computexfrstns(rmp::MasterProblem, gridtype::Symbol)
    activelines = findall(x->x!=0, round.(JuMP.getvalue(rmp.x), 5))
    linelist = rmp.linelist
    nstns = rmp.np.nstations
    stnlines = [findall(x->x!=0, [in(u, line) for line in linelist])
                for u in 1:nstns]
    xfrstops_uw = Dict{Tuple{Int,Int},Vector{Vector{Int}}}()
    xfrstops_wv = Dict{Tuple{Int,Int},Vector{Vector{Int}}}()
    for u in 1:nstns, v in nonzerodests(rmp.np, u)
        # Separate stops with active [1] and inactive [2] lines
        xfrstops_uw[u,v] = [Int[],Int[]]
        xfrstops_wv[u,v] = [Int[],Int[]]
        # Ignore all (u,v) pairs that have an active direct route.
        length(intersect(stnlines[u], stnlines[v],activelines)) > 0 && continue
        # Deal with all (u,v) pairs with a single transfer station w.
        for w in 1:nstns
            if validtransfer(rmp.np, u, v, w, rmp.angleparam, rmp.distparam, gridtype)
                uwlines = intersect(stnlines[u], stnlines[w])
                wvlines = intersect(stnlines[w], stnlines[v])
                # Lines connecting u and w are active
                if (length(uwlines) > 0) && (length(intersect(uwlines, activelines)) > 0)
                    push!(xfrstops_uw[u,v][1], w)
                # Lines connecting u and w are inactive
                # Can afford to expand set of xfrstops if nstns is small, but gets very
                # dense if nstns is large
                elseif (length(uwlines) > 0) || nstns <= 20 
                    push!(xfrstops_uw[u,v][2], w)
                end
                # Lines connecting w and v are active    
                if (length(wvlines) > 0) && (length(intersect(wvlines, activelines)) > 0)
                    push!(xfrstops_wv[u,v][1], w)
                # Lines connecting w and v are inactive 
                # Can afford to expand set of xfrstops if nstns is small, but gets very
                # dense if nstns is large
                elseif (length(wvlines) > 0) || nstns <= 20 
                    push!(xfrstops_wv[u,v][2], w)
                end
            end 
        end
    end
    xfrstops_uw, xfrstops_wv
end

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
        source_val = findfirst(round.(JuMP.getvalue(sp.src)))
        sink_val = findfirst(round.(JuMP.getvalue(sp.snk)))
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
