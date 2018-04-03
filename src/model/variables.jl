function designvariable(dm::DesignModel)
    JuMP.@variable(dm.model,
        y[1:dm.np.nsegments], Bin)
    y
end

function ridershipvariable(dm::DesignModel)
    JuMP.@variable(dm.model, 
        ride[u=1:dm.np.nstations, v=TransitNetworks.dests(dm.np,u)], Bin)
    ride 
end

function openroutevariable(dm::DesignModel)
    validod = vec([isassigned(dm.np.stage, u,v) 
        for u in 1:dm.np.nstations, v in 1:dm.np.nstations])
    stages = unique(vcat(dm.np.stage[validod]...))
    JuMP.@variable(dm.model,
        open[stages], Bin)
    open
end 
