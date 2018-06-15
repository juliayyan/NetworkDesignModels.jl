mutable struct SubProblem
    np::TN.TransitNetworkProblem
    model::JuMP.Model
    src::Vector{JuMP.Variable}
    snk::Vector{JuMP.Variable}
    edg::JuMP.JuMPDict{JuMP.Variable}
    srv::JuMP.JuMPDict{JuMP.Variable}
    srv_uw
    srv_wv
    xfrstops_uw
    xfrstops_wv
    dists::Dict{Tuple{Int,Int},Float64}
    outneighbors::Vector{Vector{Int}}
    inneighbors::Vector{Vector{Int}}
    nlegs
end

function SubProblem(
    np::TN.TransitNetworkProblem;
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    linelist::Vector{Vector{Int}} = uniquelines(np.lines),
    maxdist::Float64 = 0.5,
    direction::Vector{Float64} = [0.0,1.0],
    delta::Float64 = 1.0,
    maxlength::Int = 30,
    nlegs::Int = 1,
    transferparam::Float64 = 0.5
    )
    
    const nstns = np.nstations

    # construct graph and ensure there are no cycles
    dists = Dict{Tuple{Int,Int},Float64}() 
    outneighbors = [Int[] for u in 1:nstns]
    inneighbors  = [Int[] for u in 1:nstns]
    graph = LightGraphs.DiGraph(nstns)
    for u in 1:nstns, v in (u+1):nstns 
        d = TN.haversinedistance(np,u,v)
        b = dir(np,u,v)
        sim = dot(direction,b)/norm(direction)/norm(b)
        if (d < maxdist) && (1-abs(sim) < delta)
            if sim >= 0
                dists[u,v] = d 
                push!(outneighbors[u], v)
                push!(inneighbors[v], u)
                LightGraphs.add_edge!(graph, (u,v))
            else
                dists[v,u] = d
                push!(outneighbors[v], u)
                push!(inneighbors[u], v)    
                LightGraphs.add_edge!(graph, (v,u))
            end
        end 
    end 
    @assert !LightGraphs.is_cyclic(graph)

    # compute potential transfer stations
    if nlegs == 2
        stnlines = [find(in(u, line) for line in linelist) for u in 1:np.nstations]
        xfrstops_uw = Dict{Tuple{Int,Int},Vector{Int}}()
        xfrstops_wv = Dict{Tuple{Int,Int},Vector{Int}}()
        for u in 1:np.nstations, v in nonzerodests(np,u)
            xfrstops_uw[u,v] = Int[]
            xfrstops_wv[u,v] = Int[]
            length(intersect(stnlines[u], stnlines[v])) > 0 && continue
            for w in 1:np.nstations 
                if validtransfer(np,u,v,w,transferparam)
                    if length(intersect(stnlines[u], stnlines[w])) > 0
                        push!(xfrstops_uw[u,v], w)
                    end
                    if length(intersect(stnlines[v], stnlines[w])) > 0
                        push!(xfrstops_wv[u,v], w)
                    end
                end 
            end
        end
    else
        xfrstops_uw = xfrstops_wv = nothing
    end

    # build model
    sp = JuMP.Model(solver=solver)

    JuMP.@variable(sp, src[u=1:nstns], Bin)
    JuMP.@variable(sp, snk[u=1:nstns], Bin)
    JuMP.@variable(sp, edg[u=1:nstns, v=outneighbors[u]], Bin)
    JuMP.@variable(sp, 0 <= srv[u=1:nstns, v=nonzerodests(np,u)] <= 1)
    if nlegs == 2
        JuMP.@variable(sp, 0 <= srv_uw[u=1:nstns, v=nonzerodests(np,u)] <= 0.5)
        JuMP.@variable(sp, 0 <= srv_wv[u=1:nstns, v=nonzerodests(np,u)] <= 0.5)
        JuMP.@constraint(sp, 
            [u=1:nstns, v=nonzerodests(np,u)],
            srv[u,v] + srv_uw[u,v] + srv_wv[u,v] <= 1)
    else 
        srv_uw = srv_wv = nothing
    end

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
    if nlegs == 2
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
    end

    SubProblem(np, 
        sp, src, snk, edg, srv, srv_uw, srv_wv,
        xfrstops_uw, xfrstops_wv,
        dists, outneighbors, inneighbors,
        nlegs)
end 

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
    
    path = Int[]
    if round(JuMP.getobjectivevalue(sp.model)) > 0
        push!(path, findfirst(round.(JuMP.getvalue(sp.src))))
        count = 1
        while path[count] != findfirst(round.(JuMP.getvalue(sp.snk)))
            i = path[count]
            oldcount = count
            for j in sp.outneighbors[i]
                if round(JuMP.getvalue(sp.edg[i,j]) > 0)
                    push!(path, j)
                    count = count + 1
                    break
                end 
            end 
            if oldcount == count 
                error("Potential infinite loop")
            end
        end 
    end 
    if (length(path) > 0) && (round(sum(JuMP.getvalue(sp.edg))) != length(path)-1)
        error("Dual solution is not a valid path")
    end 
    path
end 
