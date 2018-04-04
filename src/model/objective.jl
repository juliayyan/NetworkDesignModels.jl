function ridershipobjective(dm::DesignModel)
    JuMP.@objective(dm.model, Max, sum(dm.ride))
end 
