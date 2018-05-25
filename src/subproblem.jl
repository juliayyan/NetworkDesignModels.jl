mutable struct SubProblem
    np::TN.TransitNetworkProblem
    model::JuMP.Model
    src::Vector{JuMP.Variable}
    snk::Vector{JuMP.Variable}
    edg::JuMP.JuMPDict{JuMP.Variable}
    srv::JuMP.JuMPDict{JuMP.Variable}
    dists::Dict{Tuple{Int,Int},Float64}
    outneighbors::Vector{Vector{Int}}
    inneighbors::Vector{Vector{Int}}
end

function dir(np::TN.TransitNetworkProblem, u::Int, v::Int)
    latu = np.latlngs[u,1]
    latv = np.latlngs[v,1]
    lonu = np.latlngs[u,2]
    lonv = np.latlngs[v,2]
    x = cosd(latv)*sind(lonv-lonu)
    y = cosd(latu)*sind(latv) - sind(latu)*cosd(latv)*cosd(lonv-lonu)
    [x,y]
end 

function SubProblem(
    np::TN.TransitNetworkProblem;
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    maxdist::Float64 = 0.5,
    direction::Vector{Float64} = [0.0,1.0],
    delta::Float64 = 1.0
    )
    
    const nstns = np.nstations

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

    sp = JuMP.Model(solver=solver)

    JuMP.@variable(sp, src[u=1:nstns], Bin)
    JuMP.@variable(sp, snk[u=1:nstns], Bin)
    JuMP.@variable(sp, edg[u=1:nstns, v=outneighbors[u]], Bin)
    JuMP.@variable(sp, srv[u=1:nstns, v=nonzerodests(np,u)], Bin)

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

    # demand service
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv[u,v] <= src[u] + sum(edg[w,u] for w in inneighbors[u]) + snk[u])
    JuMP.@constraint(sp,
        [u=1:nstns, v=nonzerodests(np,u)],
        srv[u,v] <= src[v] + sum(edg[w,v] for w in inneighbors[v]) + snk[v])

    SubProblem(np, sp, src, snk, edg, srv, dists, outneighbors, inneighbors)
end 

function generatecolumn(sp::SubProblem, p, q)
    JuMP.@objective(sp.model,
        Max,
        sum(sum(p[u,v]*sp.srv[u,v] for v in nonzerodests(sp.np,u)) for u in 1:sp.np.nstations) - 
        q*sum(sp.dists[u,v]*sp.edg[u,v] for u in 1:sp.np.nstations, v in sp.outneighbors[u])
    )  
    JuMP.solve(sp.model)
    
    path = Int[]
    if JuMP.getobjectivevalue(sp.model) > 0
        push!(path, findfirst(JuMP.getvalue(sp.src)))
        count = 1
        while path[count] != findfirst(JuMP.getvalue(sp.snk))
            i = path[count]
            for j in sp.outneighbors[i]
                if round(JuMP.getvalue(sp.edg[i,j]) > 0)
                    push!(path, j)
                    count = count + 1
                    break
                end 
            end 
        end 
    end 
    @assert sum(JuMP.getvalue(sp.edg)) == length(path)-1
    path
end 
