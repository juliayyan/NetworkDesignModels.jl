"""
Return a vector of stops `v` with positive demand from stop `u`.

It ensures that `np.odmatrix[u,v] > 0`, and that `u` itself does not appear.
"""
nonzerodests(np::TN.TransitNetworkProblem, u::Int) =
    filter!(i -> i != u, find(np.odmatrix[u,:] .> 0))

function uniquelines(linelist::Vector{Vector{Int}})
    uniquelines = Vector{Int}[]
    for line in linelist
        isunique = true
        for line2 in uniquelines
            if length(line) == length(line2) && all(sort(line) .== sort(line2))
                isunique = false
                break
            end
        end
        isunique && push!(uniquelines, copy(line))
    end
    uniquelines
end

linecost(np::TN.TransitNetworkProblem, line::Vector{Int}, gridtype::Symbol) =
    sum(edgecost(np, line[i], line[i+1], gridtype) for i in 1:(length(line)-1))

"""
Returns distance between stations `u` and `v` based on the `gridtype`.

If `gridtype` is `:latlong`, it returns the distance in miles. If gridtype is
`:euclidean`, it returns the euclidean distance in latlng space.
"""
function edgecost(
        np::TN.TransitNetworkProblem,
        u::Int,
        v::Int,
        gridtype::Symbol
    )
    if gridtype == :latlong
        return TN.haversinedistance(np, u, v)
    elseif gridtype == :euclidean
        return norm(np.latlngs[u,:] - np.latlngs[v,:])
    else
        error("Unrecognized gridtype: $gridtype. Use :latlong or :euclidean.")
    end
end

"""
Returns true if the commute `u` -> `w` -> `v` obeys `transferparam`.

Uses `gridtype` to determine the distances between `u,w` and `w,v` when
comparing the angle they made with transferparam.
"""
function validtransfer(
        np::TN.TransitNetworkProblem, 
        u::Int,
        v::Int,
        w::Int,
        transferparam::Float64,
        gridtype::Symbol
    )
    dir1 = dir(np, u, w, gridtype)
    dir2 = dir(np, w, v, gridtype)
    dot(dir1,dir2) / norm(dir1) / norm(dir2) >= transferparam
end

function dir(np::TN.TransitNetworkProblem, u::Int, v::Int, gridtype::Symbol)
    if gridtype == :latlong
        latu, lonu = np.latlngs[u,:]
        latv, lonv = np.latlngs[v,:]
        x = cosd(latv) * sind(lonv-lonu)
        y = cosd(latu) * sind(latv) - sind(latu) * cosd(latv) * cosd(lonv-lonu)
        return [x,y]
    elseif gridtype == :euclidean
        return np.latlngs[v,:] - np.latlngs[u,:]
    else
        error("Unrecognized gridtype: $gridtype. Use :latlong or :euclidean.")
    end
end 
