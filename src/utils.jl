nonzerodests(np::TN.TransitNetworkProblem, u::Int) =
    find(np.odmatrix[u,:] .> 0)

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

linecost(np::TN.TransitNetworkProblem, 
    line::Vector{Int}, 
    gridtype::Symbol) = 
    sum(edgecost(np, line[i], line[i+1], gridtype)
        for i in 1:length(line)-1)

function edgecost(np::TN.TransitNetworkProblem, 
    u::Int, v::Int, 
    gridtype::Symbol)
    if gridtype == :latlong
        return TN.haversinedistance(np,u,v)
    elseif gridtype == :euclidean
        return norm(np.latlngs[u,:] - np.latlngs[v,:])
    else
        @assert false
    end
end

function validtransfer(np::TN.TransitNetworkProblem, 
    u::Int, v::Int, w::Int,
    transferparam::Float64,
    gridtype::Symbol)
    dir1 = dir(np,u,w,gridtype)
    dir2 = dir(np,w,v,gridtype)
    return dot(dir1,dir2)/norm(dir1)/norm(dir2) >= transferparam                
end

function dir(np::TN.TransitNetworkProblem, u::Int, v::Int,gridtype::Symbol)
    if gridtype == :latlong
        latu = np.latlngs[u,1]
        latv = np.latlngs[v,1]
        lonu = np.latlngs[u,2]
        lonv = np.latlngs[v,2]
        x = cosd(latv)*sind(lonv-lonu)
        y = cosd(latu)*sind(latv) - sind(latu)*cosd(latv)*cosd(lonv-lonu)
        [x,y]
    elseif gridtype == :euclidean
        return np.latlngs[v,:] - np.latlngs[u,:]
    else
        @assert false
    end
end 
