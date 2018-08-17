function basemodel(
    np::TN.TransitNetworkProblem,
    inneighbors::Vector{Vector{Int}},
    outneighbors::Vector{Vector{Int}},
    maxlength::Int,
    solver)

    const nstns = np.nstations
    
    # build model
    sp = JuMP.Model(solver=solver)

    JuMP.@variable(sp, src[u=1:nstns], Bin)
    JuMP.@variable(sp, snk[u=1:nstns], Bin)
    JuMP.@variable(sp, edg[u=1:nstns, v=outneighbors[u]], Bin)
    JuMP.@variable(sp, 0 <= srv[u=1:nstns, v=nonzerodests(np,u)] <= 1)
    
    # path constraints
    JuMP.@constraint(sp,
        [u=1:nstns],
        src[u] + sum(edg[v,u] for v in inneighbors[u]) <= 1)
    JuMP.@constraint(sp,
        [u=1:nstns],
        snk[u] + sum(edg[u,v] for v in outneighbors[u]) <= 1)
    JuMP.@constraint(sp,
        [u=1:nstns],
        src[u] + sum(edg[v,u] for v in inneighbors[u]) ==
        snk[u] + sum(edg[u,v] for v in outneighbors[u]))
    JuMP.@constraint(sp,
        sum(src) == 1)
    JuMP.@constraint(sp,
        sum(snk) == 1)
    JuMP.@constraint(sp, sum(edg) <= maxlength)

    # demand service
    JuMP.@expression(sp,
        ingraph[u=1:np.nstations],
        src[u] + sum(edg[u2,u] for u2 in inneighbors[u]) + snk[u])
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv[u,v] <= ingraph[u])
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv[u,v] <= ingraph[v])

    return sp, src, snk, edg, srv, ingraph

end

function transfermodel(
    np::TN.TransitNetworkProblem,
    sp, srv, ingraph, xfrstops_uw, xfrstops_wv)

    const nstns = np.nstations

    JuMP.@variable(sp, 0 <= srv_uw[u=1:nstns, v=nonzerodests(np,u)] <= 0.5)
    JuMP.@variable(sp, 0 <= srv_wv[u=1:nstns, v=nonzerodests(np,u)] <= 0.5)
 
    JuMP.@constraint(sp, 
        [u=1:nstns, v=nonzerodests(np,u)],
        srv[u,v] + srv_uw[u,v] + srv_wv[u,v] <= 1)
    JuMP.@constraint(sp, 
        [u=1:nstns, v=nonzerodests(np,u)],
        srv[u,v] + srv_uw[u,v] + srv_wv[u,v] <= 1)
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv_uw[u,v] <= ingraph[v])
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv_wv[u,v] <= ingraph[u])
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv_uw[u,v] <= 
        ((length(xfrstops_uw[u,v]) == 0) ?
            0 : sum(ingraph[w] for w in xfrstops_uw[u,v])
            )
        )
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv_wv[u,v] <= 
        ((length(xfrstops_wv[u,v]) == 0) ?
            0 : sum(ingraph[w] for w in xfrstops_wv[u,v])
            )
        )
    
    return srv_uw, srv_wv
end