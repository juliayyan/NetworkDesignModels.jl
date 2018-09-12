"""
Constructs a base direct-route subproblem.

Use transfermodel(...) to the additional variables and constraints needed for 
modelling single-transfers.
"""
function basemodel(
        np::TN.TransitNetworkProblem,
        inneighbors::Vector{Vector{Int}},
        outneighbors::Vector{Vector{Int}},
        maxlength::Int, # maximum number of edges in a path
        solver
    )
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
        ingraph[u=1:nstns],
        src[u] + sum(edg[u2,u] for u2 in inneighbors[u]) + snk[u])
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv[u,v] <= ingraph[u])
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv[u,v] <= ingraph[v])

    return sp, src, snk, edg, srv, ingraph

end

"""
Adds transferring variables and constraints to subproblem model `sp`.

Use basemodel(...) for constructing a subproblem model.
"""
function transfermodel(
        np::TN.TransitNetworkProblem,
        sp::JuMP.Model,
        srv,
        ingraph,
        xfrstops_uw,
        xfrstops_wv
    )
    const nstns = np.nstations

    JuMP.@variable(sp, 0 <= srv_uw[u=1:nstns, v=nonzerodests(np,u)] <= 1)
    JuMP.@variable(sp, 0 <= srv_wv[u=1:nstns, v=nonzerodests(np,u)] <= 1)
 
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
            0 : sum(ingraph[w] for w in xfrstops_uw[u,v])))
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv_wv[u,v] <= 
        ((length(xfrstops_wv[u,v]) == 0) ?
            0 : sum(ingraph[w] for w in xfrstops_wv[u,v])))
    
    return srv_uw, srv_wv
end

"""
Returns a tuple (coeffs_uw, coeffs_wv) corresponding to coefficients for p.


"""
function spcoeffs(
        rmp::MasterProblem,
        sp::SubProblem
    )
    stnlines = [
        find(in(u,line) for line in rmp.linelist)
        for u in 1:rmp.np.nstations
    ]
    xval = JuMP.getvalue(rmp.x)
    coeffs_uw = Dict{Tuple{Int,Int},Float64}()
    coeffs_wv = Dict{Tuple{Int,Int},Float64}()
    for u in 1:rmp.np.nstations, v in nonzerodests(rmp.np,u)
        # all lines containing some w and either u or v
        wlines_uw = unique(vcat(union(stnlines[sp.xfrstops_uw[u,v]]...)...))
        wlines_wv = unique(vcat(union(stnlines[sp.xfrstops_wv[u,v]]...)...))
        wlines_uw = intersect(wlines_uw, stnlines[u])
        wlines_wv = intersect(wlines_wv, stnlines[v])
        if length(wlines_uw) == 0
            if length(sp.xfrstops_uw[u,v]) > 0 
                error("should have some intersection ", u, " ", v)
            end
            coeffs_uw[u,v] = 0.0
        elseif round(sum(xval[wlines_uw]),5) == 0
            coeffs_uw[u,v] = 0.5
        else 
            coeffs_uw[u,v] = 1.0
        end
        if length(wlines_wv) == 0
            if length(sp.xfrstops_wv[u,v]) > 0
                error("should have some intersection ", u, " ", v)
            end
            coeffs_wv[u,v] = 0.0
        elseif round(sum(xval[wlines_wv]),5) == 0
            coeffs_wv[u,v] = 0.5
        else 
            coeffs_wv[u,v] = 1.0
        end
    end
    coeffs_uw, coeffs_wv
end

"""
Uses dual values `p` and `q` to generate a profitable line.

This function helps solve the subproblem many times just changing the objective
without rebuilding the model.

### Keyword Arguments
* `coeffs`: a tuple (coeffs_uw, coeffs_wv) for the p variables.
* `trackingstatuses`: for tracking information through solver callbacks.

### Returns
`path`
"""
function generatecolumn(
        sp::SubProblem, 
        p,
        q;
        trackingstatuses::Vector{Symbol} = Symbol[],
        coeffs::NTuple{2, Dict{Tuple{Int,Int},Float64}} = (
            Dict(k => 0.5 for k in keys(p)),
            Dict(k => 0.5 for k in keys(p))
        )
    )
    JuMP.@objective(sp.model,
        Max,
        sum(sum(p[u,v] * (sp.srv[u,v] + (sp.nlegs == 1 ? 0 :  
                                         (coeffs[1][u,v]*sp.srv_uw[u,v] + 
                                          coeffs[2][u,v]*sp.srv_wv[u,v])))
            for v in nonzerodests(sp.np,u))
        for u in 1:sp.np.nstations) - 
        q * sum(sp.dists[u,v]*sp.edg[u,v]
            for u in 1:sp.np.nstations, v in sp.outneighbors[u])
    )

    t0 = time()
    for tracking in trackingstatuses
        if tracking == :Intermediate
            sp.auxinfo[:obj]       = Float64[]
            sp.auxinfo[:bestbound] = Float64[]
            sp.auxinfo[:time]      = Float64[]
            function boundscallback(cb)
                currobj = MathProgBase.cbgetobj(cb)
                if (length(sp.auxinfo[:obj]) == 0) || (currobj > sp.auxinfo[:obj][end])
                    push!(sp.auxinfo[:obj],       MathProgBase.cbgetobj(cb))
                    push!(sp.auxinfo[:bestbound], MathProgBase.cbgetbestbound(cb))
                    push!(sp.auxinfo[:time],      time() - t0)                
                end
            end
            JuMP.addinfocallback(sp.model, boundscallback, when = :Intermediate)
        elseif tracking == :MIPSol
            sp.auxinfo[:time] = Float64[]
            function solncallback(cb)
                push!(sp.auxinfo[:time], time() - t0)
            end
            JuMP.addinfocallback(sp.model, solncallback, when = :MIPSol)
        end
    end

    JuMP.solve(sp.model)
    sp.auxinfo[:endtime] = time() - t0

    path = getpath(sp) 
    if (length(path) > 0) && (round(sum(JuMP.getvalue(sp.edg))) != length(path)-1)
        error("Dual solution is not a valid path")
    end 

    return path
end 

"""
Heuristically generates columns.

It starts with edge-restricted networks and iteratively builds up
relevant sections of the network.

### Keyword Arguments
* `coeffs`: a tuple (coeffs_uw, coeffs_wv) for the p variables.
* `directions`: a vector of directions z to run the subproblem with
    preprocessed edge set E(z, delta).
* `nlegs`: the (maximum) number of legs of each commute.
* `maxdist`: the threshold distance for two stations to be considered
    neighbors. The distance is measured based on the `gridtype` in
    `rmp::MasterProblem`.
* `delta`: tolerance within `direction` for edges in the graph to obey.
* `maxlength`: the maximum number of edges in a path
* `stepsize`:
* `solver`: The solver being used to solve the problem.
* `maxiterations`: The maximum number of iterations to build up relevant
        sections of the network for.
* `trackingstatuses`: for tracking information through solver callbacks.

### Returns
A `(path, auxinfo)` tuple.
"""
function generatecolumn(
        rmp::MasterProblem; 
        coeffs::NTuple{2, Dict{Tuple{Int,Int},Float64}} = (
            Dict(k => 0.5 for k in keys(p)),
            Dict(k => 0.5 for k in keys(p))
        ),
        directions::Vector{Vector{Float64}} = [
            [ 1.0, 0.0], [ 1.0, 0.5],
            [ 1.0, 1.0], [ 0.5, 1.0],
            [ 0.0, 1.0], [-0.5, 1.0],
            [-1.0, 1.0], [-1.0,-0.5]
        ],
        nlegs::Int = length(rmp.commutelines),
        maxdist::Float64 = 0.5,
        delta::Float64 = 1.0,
        maxlength::Int = 30,
        solver = Gurobi.GurobiSolver(OutputFlag = 0),
        stepsize::Int = 1,
        maxiterations::Int = 10,
        trackingstatuses::Vector{Symbol} = Symbol[]
    )
    auxinfo = Dict{Symbol,Any}()
    auxinfo[:time]   = Float64[]
    auxinfo[:obj]    = Float64[]
    auxinfo[:nnodes] = Int[]
    auxinfo[:nlazy]  = 0
    t0 = time()

    const nstns = rmp.np.nstations
    const dists = [edgecost(rmp.np,u,v,rmp.gridtype) for u in 1:nstns, v in 1:nstns]
    const neighbors = [setdiff(find(dists[u,:] .< maxdist),u) for u in 1:nstns]

    p = JuMP.getdual(rmp.choseline)
    q = JuMP.getdual(rmp.bcon)
    # s = JuMP.getdual(rmp.choseub)

    # compute warm start
    oldobj = 0.0
    sp_warm = [SubProblem(rmp,
                          nlegs = nlegs,
                          maxdist = maxdist, 
                          maxlength = maxlength,
                          direction = d, 
                          solver = Gurobi.GurobiSolver(OutputFlag = 0)) 
                for d in directions]
    sp_objs = [JuMP.getobjectivevalue(spw.model) for spw in sp_warm]
    soln_warm = [generatecolumn(spw,p,q,coeffs=coeffs) for spw in sp_warm]
    nodeset = unique(union(soln_warm...))
    push!(auxinfo[:time], time() - t0)
    push!(auxinfo[:obj], maximum(sp_objs))
    push!(auxinfo[:nnodes], length(nodeset))

    # iteratively generate path
    path = soln_warm[findmax(sp_objs)[2]]
    for i in 1:maxiterations
        sp = SubProblemCP(rmp, 
                          nlegs = nlegs,
                          maxdist = maxdist, 
                          maxlength = maxlength,
                          nodeset = nodeset, 
                          solver = solver)
        warmstart(sp, path)
        path = generatecolumn(sp,
                              p,
                              q,
                              coeffs=coeffs,
                              trackingstatuses=trackingstatuses)
        if round(JuMP.getobjectivevalue(sp.model),3) == round(oldobj,3)
            break
        else
            auxinfo[:nlazy] += sp.auxinfo[:nlazy]
            push!(auxinfo[:time], time() - t0)
            push!(auxinfo[:obj],  JuMP.getobjectivevalue(sp.model))
            nodeset = union(nodeset, path)
            for j in 2:stepsize 
                nodeset = union(nodeset, neighbors[nodeset]...)
            end
            push!(auxinfo[:nnodes], length(nodeset))
            oldobj = JuMP.getobjectivevalue(sp.model)
        end
    end

    auxinfo[:endtime] = time() - t0
    return path, auxinfo 
end

function warmstart(sp::SubProblem, path::Vector{Int})
    JuMP.setvalue(sp.src[path[1]]  , 1)
    JuMP.setvalue(sp.snk[path[end]], 1)
    for i in 2:length(path)
        JuMP.setvalue(sp.edg[path[i-1],path[i]], 1)
    end
end
