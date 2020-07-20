# SECTION 4.2.2 BLUE BIKES (formerly Hubway)

using NetworkDesignModels, Gurobi, JuMP
using DataFrames, JLD2, CSV
using ProgressMeter
NDM = NetworkDesignModels;

# parameters
nzs = 4

# load data
env = Gurobi.Env()
options = NDM.MasterOptions(
    nlegs = 2, distparam = 1.5,
    nfreqs = 2, costwts = [1.0, 0.3], freqwts = [1.0, 0.45], xfrwts = [1.0, 0.0], 
    constrainedg = true, ninitial = 9)
JLD2.@load "data/transit-network-boston-hubway.jld2" stns np
@show sum(np.odmatrix .> 0)

for budget in [30, 40, 50]

# create master problem
options.modeltype = :lp
rmp = NDM.MasterProblem(
    np, 
    linelist = copy(np.lines),
    options = options,
    solver = GurobiSolver(env, OutputFlag = 0)
);
NDM.optimize(rmp, budget)
@show length(NDM.commutes(rmp))

directions = [round.([cos(f*pi),sin(f*pi)], digits=3) for f in (1/nzs):(1/nzs):1]
curtime = time()
for it in 1:20

    @show it
    NDM.optimize(rmp, budget)
    println(getobjectivevalue(rmp.model))
    pathstoadd = Vector{Vector{Int}}()

    for f in 1:options.nfreqs

        objs = Vector{Float64}()
        paths = Vector{Vector{Int}}()
        lastpath = nothing
        # solve
        @showprogress for d in 1:length(directions)
            sp_d = NDM.SubProblem(rmp,
                maxdist = 1.0,
                solver = GurobiSolver(env, OutputFlag = 0, TimeLimit = 1*60),
                maxlength = 19,
                direction = directions[d], delta = 1.0, # 1 - (delta in paper) -- 1.0 in the code corresponds to 0.0 in the paper
                traveltimes = true);
            @constraint(sp_d.model, sum(sp_d.edg) >= 7)
            #=if lastpath != nothing
                NDM.trywarmstart(sp_d, lastpath)
            end=#
            path = NDM.generatecolumn(sp_d, rmp, f=f)
            if length(path) > 0 && !NDM.linein(path, rmp.linelist)
                push!(paths, path)
                push!(objs, getobjectivevalue(sp_d.model))
                lastpath = path
            end
        end
        if length(objs) == 0
            continue
        else
            path_warm = paths[findmax(objs)[2]]
            push!(pathstoadd, path_warm)
        end

    end

    added = 0
    if length(pathstoadd) > 0
        for path_warm in pathstoadd
            if !NDM.linein(path_warm, rmp.linelist)
                NDM.addcolumn!(rmp, path_warm)
                added += 1
            end
        end
    end
    if added == 0
        break
    end

end
runtime = time() - curtime

# solve IP
options.modeltype = :ip
rmpip0 = NDM.MasterProblem(
    np, 
    linelist = copy(np.lines),
    options = options,
    solver = GurobiSolver(env, OutputFlag = 0)
);
NDM.optimize(rmpip0, budget)
println(getobjectivevalue(rmpip0.model))
rmpip = NDM.MasterProblem(
    np, 
    linelist = rmp.linelist,
    options = options,
    solver = GurobiSolver(env, OutputFlag = 0)
);
NDM.optimize(rmpip, budget)
println(getobjectivevalue(rmpip.model))

twolines = 0
total = 0
for (u,v) in NDM.commutes(np)
    uv1 = getvalue(rmpip.θ[(u,v)])
    if length(np.xfrstns[(u,v)]) > 0
        uv2 = getvalue(sum(rmpip.model[:twoline][(u,v),w] for w in np.xfrstns[(u,v)]))
        twolines += uv2
    end
    total += uv1
end
println(twolines/total)

df = DataFrame(line = Int[], order = Int[], station = Int[], lat = Float64[], lon = Float64[], freq = Int[], old = Int[])
for l in 1:length(rmp.linelist)
    for f in 1:options.nfreqs
        if getvalue(rmpip.x[l,f]) > 0.1
            line = rmp.linelist[l]
            for i in 1:length(line)
                u = line[i]
                push!(df, [l, i, u, np.latlon[u, 1], np.latlon[u, 2], f, 1*(l <= length(np.lines))])
            end
        end
    end
end

# JLD2.@load "data/boston-hubway-all-results.jld2" allbudgets alloptions allcolumns allruntimes
# push!(allbudgets, budget)
# push!(alloptions, options)
# push!(allcolumns, rmp.linelist)
# push!(allruntimes, runtime)
# JLD2.@save "data/boston-hubway-all-results.jld2" allbudgets alloptions allcolumns allruntimes

end