"""
SubProblem
"""
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

"""
SubProblem Constructor
* Constructs an acyclic graph that only has edges 
* within some `delta`tolerance of `direction`
"""
function SubProblem(
    rmp::MasterProblem;
    nlegs::Int = length(rmp.commutelines), # want this in case 2-legs hard to solve
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    maxdist::Float64 = 0.5,
    direction::Vector{Float64} = [0.0,1.0],
    delta::Float64 = 1.0,
    maxlength::Int = 30
    )
    
    const np = rmp.np
    const nstns = np.nstations
    const gridtype = rmp.gridtype

    # construct graph and ensure there are no cycles
    dists = Dict{Tuple{Int,Int},Float64}() 
    outneighbors = [Int[] for u in 1:nstns]
    inneighbors  = [Int[] for u in 1:nstns]
    graph = LightGraphs.DiGraph(nstns)
    for u in 1:nstns, v in (u+1):nstns 
        d = edgecost(np, u, v, gridtype)
        b = dir(np,u,v,gridtype)
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

    sp, src, snk, edg, srv, ingraph = basemodel(np, inneighbors, outneighbors, maxlength, solver)

    if nlegs == 2
        xfrstops_uw, xfrstops_wv = computexfrstns(np, rmp.linelist, rmp.transferparam, gridtype)
        srv_uw, srv_wv = transfermodel(np, sp, srv, ingraph, xfrstops_uw, xfrstops_wv)
    else 
        xfrstops_uw = xfrstops_wv = nothing
        srv_uw = srv_wv = nothing
    end

    SubProblem(np, 
        sp, src, snk, edg, srv, srv_uw, srv_wv,
        xfrstops_uw, xfrstops_wv,
        dists, outneighbors, inneighbors,
        nlegs)
end 

"uses dual values `p`,`q`,`s` to generate a profitable line"
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
