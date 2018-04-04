function routeopenconstraint(dm::DesignModel)
    validod = vec([isassigned(dm.np.stage, u,v) 
        for u in 1:dm.np.nstations, v in 1:dm.np.nstations])
    stages = unique(vcat(dm.np.stage[validod]...))
    JuMP.@constraint(dm.model, 
        [r=1:length(stages), s=stages[r]], 
        dm.open[stages[r]] <= dm.y[s])    
    JuMP.@constraint(dm.model,
        [r=1:length(stages)],
        dm.open[stages[r]] >=
        sum(dm.y[s] for s in stages[r]) - length(stages[r]) + 1)
end

function ridetransitconstraint(dm::DesignModel)
    routeoptions = Dict{Tuple{Int, Int}, Any}()
    for u=1:dm.np.nstations, v=TransitNetworks.dests(dm.np, u)
        # TODO: Replace with stages that are ranked higher than 
        # no-transit option in preference ranking
        routeoptions[u,v] = 
            unique(TransitNetworks.expandedstages(dm.np,u,v)[2])
    end 
    JuMP.@constraint(dm.model,
        [u=1:dm.np.nstations, 
         v=TransitNetworks.dests(dm.np,u),
         r=1:length(routeoptions[u,v])],
        dm.ride[u,v] >= dm.open[routeoptions[u,v][r]])
    JuMP.@constraint(dm.model,
        [u=1:dm.np.nstations,
         v=TransitNetworks.dests(dm.np,u)],
        dm.ride[u,v] <= 
        sum(dm.open[routeoptions[u,v][r]] 
            for r in 1:length(routeoptions[u,v])))
end 

function budgetconstraint(dm::DesignModel)
    costs = [sum([TransitNetworks.haversinedistance(dm.np, 
                    seg[i], seg[i+1]) for i in 1:length(seg)-1])
                for seg in dm.np.segments]
    JuMP.@constraint(dm.model,
        dot(costs, dm.y) <= dm.np.budget)
end 
