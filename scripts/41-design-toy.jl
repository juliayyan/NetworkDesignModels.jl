# SECTION 4.1

using NetworkDesignModels, Gurobi, JuMP
using DataFrames, JLD2

# load data
env = Gurobi.Env()
results = DataFrame(
    budget = Float64[], niterations = Int[], runtime = Float64[],
    line = Int[], order = Int[], stop_id = String[], freq = Float64[], 
    pctxfr = Float64[], obj = Float64[])
JLD2.@load "data/transit-network-toy.jld2" stns np

# create master problem
options = NetworkDesignModels.MasterOptions(
    nlegs = 2, 
    nfreqs = 2, costwts = [1.0, 0.5], freqwts = [1.0, 0.75], xfrwts = [1.0, 0.25])
rmp = NetworkDesignModels.MasterProblem(
    np, 
    linelist = Vector{Int}[],
    options = options,
    solver = GurobiSolver(env, OutputFlag = 0)
);
scale = 1.0
budget = ceil(scale*sum(NetworkDesignModels.linecost(np,line) for line in np.lines))
NetworkDesignModels.optimize(rmp, budget)

# create subproblem objects
sp1 = NetworkDesignModels.SubProblemCP(rmp, maxdist = 1.45, traveltimes = true, 
    solver = GurobiSolver(env, OutputFlag = 0))
sp2 = NetworkDesignModels.SubProblemCP(rmp, nlegs = 1, maxdist = 1.45, traveltimes = true, 
    solver = GurobiSolver(env, OutputFlag = 0))

# generate columns
t0 = time()
nit = 0
for i in 1:30

    global nit += 1

    # solve master
    NetworkDesignModels.optimize(rmp, budget)
    println(getobjectivevalue(rmp.model))
    
    # create subproblem
    global sp1 = NetworkDesignModels.SubProblemCP(rmp, maxdist = 1.45, traveltimes = true, 
        solver = GurobiSolver(env, OutputFlag = 0))
    path = NetworkDesignModels.generatecolumn(sp1, rmp)
    
    if length(path) > 0 && !NetworkDesignModels.linein(path, rmp.linelist)
        println("$(i) adding high frequency column")
        for k1 in 1:(length(path)-1), k2 in (k1+1):length(path)
            if !NetworkDesignModels.linein(path[k1:k2], rmp.linelist)
                NetworkDesignModels.addcolumn!(rmp, path[k1:k2])
            end
        end
    elseif length(path) > 0
        error("Bad column")
        break
    else
        global sp2 = NetworkDesignModels.SubProblemCP(rmp, nlegs = 1, maxdist = 1.45, traveltimes = true, 
            solver = GurobiSolver(env, OutputFlag = 0))
        path = NetworkDesignModels.generatecolumn(sp2, rmp, f=2)
        if length(path) > 0 && !NetworkDesignModels.linein(path, rmp.linelist)
            println("$(i) adding low frequency column")
            for k1 in 1:(length(path)-1), k2 in (k1+1):length(path)
                if !NetworkDesignModels.linein(path[k1:k2], rmp.linelist)
                    NetworkDesignModels.addcolumn!(rmp, path[k1:k2])
                end
            end
        elseif length(path) > 0
            error("bad column ", path)
            break
        else
            break
        end
    end
end
t0 = time() - t0

# solve IP
options.modeltype = :ip
rmpip = NetworkDesignModels.MasterProblem(
    np, 
    linelist = rmp.linelist,
    options = options,
    solver = GurobiSolver(env, OutputFlag = 0)
);
NetworkDesignModels.optimize(rmpip, budget)
println(getobjectivevalue(rmpip.model))

# save data
twolines = 0
total = 0
for (u,v) in NetworkDesignModels.commutes(np)
    uv1 = getvalue(rmpip.θ[(u,v)])
    if length(np.xfrstns[(u,v)]) > 0
        uv2 = getvalue(sum(rmpip.model[:twoline][(u,v),w] for w in np.xfrstns[(u,v)]))
        global twolines += uv2
    end
    global total += uv1
end
println(twolines/total)

finallist = Vector{Vector{Int}}()
for l in 1:length(rmpip.linelist), f in 1:options.nfreqs
    if getvalue(rmpip.x[l,f]) > 1e-3
        println(f, "\t", rmpip.linelist[l])
        push!(finallist, rmpip.linelist[l])
        for j in 1:length(rmpip.linelist[l])
            push!(results, [budget, nit, t0,
                l, j, stns[rmpip.linelist[l][j]], f,
                twolines/total,
                getobjectivevalue(rmpip.model)])
        end
    end
end

