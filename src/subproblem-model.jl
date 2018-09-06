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

    JuMP.@variable(sp, 0 <= srv2[u=1:nstns, v=nonzerodests(np,u)] <= 1)
 
    JuMP.@constraint(sp, 
        [u=1:nstns, v=nonzerodests(np,u)],
        srv[u,v] + srv2[u,v] <= 1)
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv2[u,v] <= ingraph[u] + ingraph[v])
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv2[u,v] <= 
        ((length(xfrstops_uw[u,v]) == 0) ?
            0 : sum(ingraph[w] for w in union(xfrstops_uw[u,v],xfrstops_wv[u,v]))
            )
        )
    
    return srv2
end

"uses dual values `p`,`q`,`s` to generate a profitable line.
 this function helps solve the subproblem many times
 without rebuilding the model, just changing the objective.
"
function generatecolumn(sp::SubProblem, 
    p, q, s;
    tracking::Symbol = :none)
    JuMP.@objective(sp.model,
        Max,
        sum(sum(min(p[u,v],sp.np.odmatrix[u,v] - s[u,v])*
            (sp.srv[u,v] + 
                (sp.nlegs == 1 ? 0 : 
                 0.5*sp.srv2[u,v])) 
            for v in nonzerodests(sp.np,u)
                )
        for u in 1:sp.np.nstations) - 
        q*sum(sp.dists[u,v]*sp.edg[u,v] 
            for u in 1:sp.np.nstations, 
                v in sp.outneighbors[u])
    )  

    t0 = time()

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

    JuMP.solve(sp.model)
    sp.auxinfo[:endtime] = time() - t0

    path = getpath(sp) 
    if (length(path) > 0) && (round(sum(JuMP.getvalue(sp.edg))) != length(path)-1)
        error("Dual solution is not a valid path")
    end 

    return path
end 

"Heuristically generates columns by starting with edge-restricted
 networks and iteratively building up relevant sections of network"
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
    stepsize::Int = 1,
    maxiterations::Int = 10,
    tracking::Symbol = :none
    )
    
    warn("Only implemented for nlegs = 1")

    auxinfo = Dict{Symbol,Any}()
    auxinfo[:time]   = Float64[]
    auxinfo[:obj]    = Float64[]
    auxinfo[:nnodes] = Int[]
    auxinfo[:nlazy]  = 0
    t0 = time()

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
    push!(auxinfo[:time], time() - t0)
    push!(auxinfo[:obj], maximum(sp_objs))
    oldobj = 0
    nodeset = soln_warm[findmax(sp_objs)[2]]
    push!(auxinfo[:nnodes], length(nodeset))

    # iteratively generate path
    path = nodeset
    for i = 1:maxiterations
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
