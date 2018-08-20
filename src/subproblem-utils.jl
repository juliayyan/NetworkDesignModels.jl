"computes all potential transfer stations between commutes (u,v)"
function computexfrstns(
    np::TN.TransitNetworkProblem,
    linelist::Vector{Vector{Int}},
    transferparam::Float64,
    gridtype::Symbol
)
    stnlines = [find(in(u, line) for line in linelist) for u in 1:np.nstations]
    xfrstops_uw = Dict{Tuple{Int,Int},Vector{Int}}()
    xfrstops_wv = Dict{Tuple{Int,Int},Vector{Int}}()
    for u in 1:np.nstations, v in nonzerodests(np,u)
        xfrstops_uw[u,v] = Int[]
        xfrstops_wv[u,v] = Int[]
        length(intersect(stnlines[u], stnlines[v])) > 0 && continue
        for w in 1:np.nstations 
            if validtransfer(np,u,v,w,transferparam,gridtype)
                if length(intersect(stnlines[u], stnlines[w])) > 0
                    push!(xfrstops_uw[u,v], w)
                end
                if length(intersect(stnlines[v], stnlines[w])) > 0
                    push!(xfrstops_wv[u,v], w)
                end
            end 
        end
    end
    return xfrstops_uw, xfrstops_wv
end

"finds path from source to sink using subproblem `sp`"
function getpath(sp::SubProblem)
    if round(JuMP.getobjectivevalue(sp.model)) > 0
        np = sp.np
        visited = falses(np.nstations)
        source_val = findfirst(round.(JuMP.getvalue(sp.src)))
        sink_val = findfirst(round.(JuMP.getvalue(sp.snk)))
        return getpath(source_val, 
            source_val, sink_val,
            sp.edg,
            visited, sp.outneighbors)
    end
    return Int[]
end

"helper function for wrapper getpath(), where `cur` can be 
 any starting node (in case of subtours)"
function getpath(cur::Int, 
    source_val::Int, sink_val::Int, 
    edg::JuMP.JuMPDict{JuMP.Variable}, 
    visited::BitArray,
    outneighbors::Vector{Vector{Int}};
    maxiterations = 1000)
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
