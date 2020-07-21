# SECTION 4.2.1

using NetworkDesignModels, Gurobi, JuMP
using DataFrames, JLD2, CSV
using ProgressMeter
NDM = NetworkDesignModels;

# parameters
fulltime = 20
nzs = 4
if length(ARGS) == 0
    delta = 1.0 # 1 - (delta in paper) -- 1.0 in the code corresponds to 0.0 in the paper
    distparam = 1.5
else
    delta = [1.0, 0.9, 0.8][parse(Int, ARGS[1])]
    distparam = [1.4, 1.5][parse(Int, ARGS[2])]
end
@show delta
@show distparam

# load data
env = Gurobi.Env()
options = NDM.MasterOptions(
    nlegs = 1, distparam = distparam,
    nfreqs = 1, costwts = [1.0], freqwts = [1.0], xfrwts = [1.0])
JLD2.@load "data/transit-network-boston.jld2" stns np
@show sum(np.odmatrix .> 0)

# create master problem
rmp = NDM.MasterProblem(
    np, 
    linelist = copy(np.lines),
    options = options,
    solver = GurobiSolver(env, OutputFlag = 0)
);
budget = ceil(0.75*sum(NDM.linecost(np,line) for line in np.lines))
NDM.optimize(rmp, budget)
@show length(NDM.commutes(rmp))
@time badcoms = NDM.badcommutecombos(np, options, k=2, rmp=rmp)
@show length(badcoms)

# generate warm start from processed subproblems
directions = [round.([cos(f*pi),sin(f*pi)], digits=3) for f in (1/nzs):(1/nzs):1]
objs = Vector{Float64}()
paths = Vector{Vector{Int}}()
lastpath = nothing
# warm up functions to get accurate running times
sp_d = NDM.SubProblem(rmp,
        maxdist = 1.0,
        solver = GurobiSolver(env, OutputFlag = 0, TimeLimit = 5),
        maxlength = 19,
        direction = directions[1], delta = delta,
        traveltimes = true)
JuMP.@constraint(sp_d.model, 
    [bc in badcoms], 
    sum(sp_d.srv[uv] for uv in bc) <= length(bc)-1)
temp = NDM.generatecolumn(sp_d, rmp)
if lastpath != nothing
    NDM.trywarmstart(sp_d, temp)
end
# solve
time_warm = time()
for d in 1:length(directions)
    global sp_d = NDM.SubProblem(rmp,
        maxdist = 1.0,
        solver = GurobiSolver(env, OutputFlag = 0, TimeLimit = 1*60),
        maxlength = 19,
        direction = directions[d], delta = delta,
        traveltimes = true);
    # cuts
    JuMP.@constraint(sp_d.model, 
        [bc in badcoms], 
        sum(sp_d.srv[uv] for uv in bc) <= length(bc)-1)
    if lastpath != nothing
        NDM.trywarmstart(sp_d, lastpath)
    end
    path = NDM.generatecolumn(sp_d, rmp)
    if length(path) > 0 && !NDM.linein(path, rmp.linelist)
        push!(paths, path)
        push!(objs, getobjectivevalue(sp_d.model))
        global lastpath = path
    end
end
time_warm = time() - time_warm
obj_warm = maximum(objs)
println(obj_warm)
path_warm = paths[findmax(objs)[2]]

if delta == 1.0

# solve full subproblem
sp = NDM.SubProblemCP(rmp,
    maxdist = 1.0,
    solver = GurobiSolver(env, OutputFlag = 1, TimeLimit = fulltime*60),
    maxlength = 19,
    traveltimes = true);
NDM.trywarmstart(sp, path_warm)
# cut off reverse solution because this graph is symmetric
for k in 2:length(path_warm)
    JuMP.@constraint(sp.model, sp.edg[path_warm[k],path_warm[k-1]] == 0)
end
# cuts
JuMP.@constraint(sp.model, 
    [bc in badcoms], 
    sum(sp.srv[uv] for uv in bc) <= length(bc)-1)
path_full = NDM.generatecolumn(sp, rmp, trackingstatuses = [:Intermediate])
time_full = sp.auxinfo[:endtime]
obj_full = getobjectivevalue(sp.model)
time_progress_full = sp.auxinfo[:time]
obj_progress_full = sp.auxinfo[:obj]
bound_progress_full = sp.auxinfo[:bestbound]
gap_full = sp.auxinfo[:bestbound][end]/obj_full-1

# save data
# JLD2.@save "data/results-scale-$(delta)-$(nzs)-$(distparam).jld2" paths objs time_warm obj_warm path_warm time_full obj_full path_full gap_full time_progress_full obj_progress_full bound_progress_full

else 

# save data
# JLD2.@save "data/results-scale-$(delta)-$(nzs)-$(distparam).jld2" paths objs time_warm obj_warm path_warm

end