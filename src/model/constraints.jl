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
