"""
Return a vector of stops `v` with positive demand from stop `u`.

It ensures that `np.odmatrix[u,v] > 0`, and that `u` itself does not appear.
"""
nonzerodests(np::TransitNetwork, u::Int) =
    filter!(i -> i != u, findall(np.odmatrix[u,:] .> 0))

"""
Return a vector of commutes with nonzero demand.
"""
commutes(np::TransitNetwork) =
    unique(
        vcat([[(min(u,v),max(u,v)) for v in nonzerodests(np,u)] for u in 1:np.nstations]...)
    )

"""
Return a vector of commutes with demand > lb
"""
function commutes(np::TransitNetwork, lb::Int)
    coms = Vector{Tuple{Int,Int}}()
    for (u,v) in commutes(np)
        if demand(np, (u,v)) > lb
            push!(coms, (u,v))
        end
    end
    coms
end

"""
Return a vector of commutes with nonzero dual variables.
"""
function commutes(rmp::MasterProblem) 
    p = JuMP.getdual(rmp.model[:choseline])
    alluvs = commutes(rmp.np)
    smalluvs = Vector{Tuple{Int,Int}}()
    for (u,v) in alluvs
        if p[(u,v)] > 0
            push!(smalluvs, (u,v))
        end
    end
    smalluvs
end

"Total demand between u and v"
demand(np::TransitNetwork, uv::Tuple{Int,Int}) =
    np.odmatrix[uv[1],uv[2]] + np.odmatrix[uv[2],uv[1]]

"""
Return the unique lines of `linelist`.

It ensures that both a line and its reverse will not appear.
"""
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

"""
Return whether `line` or its reversal appears in `linelist`.
"""
linein(line::Vector{Int}, linelist::Vector{Vector{Int}}) = in(line, linelist) || in(reverse(line), linelist)

"""
The total distance along a line
"""
linecost(np::TransitNetwork, line::Vector{Int}) =
    sum(edgecost(np, line[i], line[i+1]) for i in 1:(length(line)-1))

"""
Returns distance between stations `u` and `v` based on the `gridtype`.

If `gridtype` is `:latlong`, it returns the distance in miles. If gridtype is
`:euclidean`, it returns the euclidean distance in latlng space.
"""
function edgecost(
        np::TransitNetwork,
        u::Int,
        v::Int
    )
    if np.dists != nothing
        return haskey(np.dists, (u,v)) ? np.dists[u,v] : Inf 
    elseif np.gridtype == :latlong
        return TN.haversinedistance(np, u, v)
    elseif np.gridtype == :euclidean
        return norm(np.latlon[u,:] - np.latlon[v,:])
    else
        error("Unrecognized gridtype: $np.gridtype. Use :latlong or :euclidean.")
    end
end

"""
Returns the distance of the shortest path between u and v.
If distance is not specified in TransitNetwork, calculates a default
value based on directly going from u to v via haversine or euclidean distance.
"""
function spcost(
    np::TransitNetwork,
    u::Int,
    v::Int)
    if np.spdists != nothing
        return haskey(np.spdists, (u,v)) ? np.spdists[u,v] : Inf
    else
        return edgecost(np,u,v)
    end
end

"""
Returns the direction of edge (u, v) projected onto cartesian coordinates
"""
function dir(np::TransitNetwork, u::Int, v::Int)
    if np.gridtype == :latlong
        latu, lonu = np.latlon[u,:]
        latv, lonv = np.latlon[v,:]
        x = cosd(latv) * sind(lonv-lonu)
        y = cosd(latu) * sind(latv) - sind(latu) * cosd(latv) * cosd(lonv-lonu)
        return [x,y]
    elseif gridtype == :euclidean
        return np.latlon[v,:] - np.latlon[u,:]
    else
        error("Unrecognized gridtype: $gridtype. Use :latlong or :euclidean.")
    end
end 
