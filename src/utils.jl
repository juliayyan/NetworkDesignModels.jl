nonzerodests(np::TN.TransitNetworkProblem, u::Int) =
    find(np.odmatrix[u,:] .> 0)

linecost(np::TN.TransitNetworkProblem, line::Vector{Int}) = 
    sum(TN.haversinedistance(np, line[i], line[i+1])
        for i in 1:length(line)-1)

function uniquelines(linelist::Vector{Vector{Int}})
    uniquelines = Vector{Int}[]
    for line in linelist
        isunique = true
        for line2 in uniquelines 
            if (length(line) == length(line2)) && 
               all(sort(line) .== sort(line2))
                isunique = false
                break
            end
        end
        isunique && push!(uniquelines, line)
    end
    uniquelines
end

function validtransfer(np::TN.TransitNetworkProblem, 
    u::Int, v::Int, w::Int;
    maxdot::Float64 = 0.5)
    dir1 = dir(np,u,w)
    dir2 = dir(np,w,v)
    return dot(dir1,dir2)/norm(dir1)/norm(dir2) >= maxdot                
end