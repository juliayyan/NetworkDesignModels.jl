"base direct-route subproblem"
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
    JuMP.@constraint(sp, [u=1:nstns], src[u] + snk[u] <= 1)
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

"adds transferring variables and constraints to subproblem model `sp`"
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

"uses dual values `p`,`q`,`s` to generate a profitable line.
 this function helps solve the subproblem many times
 without rebuilding the model, just changing the objective.
"
function generatecolumn(sp::SubProblem, p, q, s)
    JuMP.@objective(sp.model,
        Max,
        sum(sum(min(p[u,v],sp.np.odmatrix[u,v] - s[u,v])*
            (sp.srv[u,v] + 
                (sp.nlegs == 1 ? 0 : 
                 sp.srv_uw[u,v] + sp.srv_wv[u,v])) 
            for v in nonzerodests(sp.np,u)
                )
        for u in 1:sp.np.nstations) - 
        q*sum(sp.dists[u,v]*sp.edg[u,v] 
            for u in 1:sp.np.nstations, 
                v in sp.outneighbors[u])
    )  
    JuMP.solve(sp.model)
    path = getpath(sp) 
    if (length(path) > 0) && (round(sum(JuMP.getvalue(sp.edg))) != length(path)-1)
        error("Dual solution is not a valid path")
    end 
    path
end 

""
function generatecolumn(rmp::MasterProblem; 
    directions::Vector{Vector{Float64}} = [
        [1.0,0.0],[1.0,0.5],
        [1.0,1.0],[0.5,1.0],
        [0.0,1.0],[-0.5,1.0],
        [-1.0,1.0],[-1.0,-0.5]
    ],
    maxdist::Float64 = 0.5,
    maxlength::Int = 30,
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    stepsize::Int = 1)

    const np = rmp.np
    const dists = [edgecost(np,u,v,rmp.gridtype) for u in 1:np.nstations, v in 1:np.nstations]
    const neighbors = [setdiff(find(dists[u,:] .< maxdist),u) for u in 1:np.nstations]

    p = JuMP.getdual(rmp.choseline)
    q = JuMP.getdual(rmp.bcon)
    s = JuMP.getdual(rmp.choseub)

    # compute warm start
    sp_warm = [SubProblem(rmp,
                          maxdist = maxdist, 
                          maxlength = maxlength,
                          direction = d, 
                          solver = Gurobi.GurobiSolver(OutputFlag = 0)) 
                for d in directions]
    soln_warm = [generatecolumn(spw,p,q,s) for spw in sp_warm]
    sp_objs = [JuMP.getobjectivevalue(spw.model) for spw in sp_warm]
    @show sp_objs
    oldobj = 0
    nodeset = soln_warm[findmax(sp_objs)[2]]

    # iteratively generate path
    path = nodeset
    for i = 1:10
        @show length(nodeset)
        sp = SubProblemCP(rmp, 
                          maxdist = maxdist, 
                          maxlength = maxlength,
                          nodeset = nodeset, 
                          solver = solver)
        warmstart(sp, path)
        path = generatecolumn(sp,p,q,s)
        if round(JuMP.getobjectivevalue(sp.model),3) == round(oldobj,3)
            break
        else
            nodeset = union(nodeset, path)
            for j in 2:stepsize 
                nodeset = union(nodeset, neighbors[nodeset]...)
                @show length(nodeset)
            end
            oldobj = JuMP.getobjectivevalue(sp.model)
        end
    end

    path
end

function warmstart(sp::SubProblem, path::Vector{Int})
    JuMP.setvalue(sp.src[path[1]]  , 1)
    JuMP.setvalue(sp.snk[path[end]], 1)
    for i in 2:length(path)
        JuMP.setvalue(sp.edg[path[i-1],path[i]], 1)
    end
end
