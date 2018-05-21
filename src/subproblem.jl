mutable struct SubProblem
    np::TN.TransitNetworkProblem
    model::JuMP.Model
    src::Vector{JuMP.Variable}
    snk::Vector{JuMP.Variable}
    edg::JuMP.JuMPDict{JuMP.Variable}
    srv::JuMP.JuMPDict{JuMP.Variable}
    dists::Matrix{Float64}
    outneighbors::Vector{Vector{Int}}
    inneighbors::Vector{Vector{Int}}
end

function SubProblem(
    np::TN.TransitNetworkProblem;
    solver = Gurobi.GurobiSolver(OutputFlag = 0),
    maxdist::Float64 = 0.5
    )
    
    const nstns = np.nstations

    dists = zeros(nstns, nstns); 
    for u in 1:nstns, v in (u+1):nstns 
        dists[u,v] = TN.haversinedistance(np,u,v)
        dists[v,u] = dists[u,v]
    end 
    outneighbors = [setdiff(find(dists[i,:] .< maxdist), 1:i) for i in 1:nstns]
    inneighbors = [Int[] for u in 1:nstns]
    for u in 1:nstns
        n = length(outneighbors[u])
        for i in 1:n 
            v = outneighbors[u][i]
            push!(inneighbors[v], u)
        end 
    end

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
    path
end 
