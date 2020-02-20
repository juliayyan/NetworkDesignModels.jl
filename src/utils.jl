"""
Return a vector of stops `v` with positive demand from stop `u`.

It ensures that `np.odmatrix[u,v] > 0`, and that `u` itself does not appear.
"""
nonzerodests(np::TN.TransitNetworkProblem, u::Int) =
    filter!(i -> i != u, findall(np.odmatrix[u,:] .> 0))

"""
Return a vector of commutes with nonzero demand.
"""
commutes(np::TN.TransitNetworkProblem) =
    unique(
        vcat([[(min(u,v),max(u,v)) for v in nonzerodests(np,u)] for u in 1:np.nstations]...)
    )

demand(np::TN.TransitNetworkProblem, uv::Tuple{Int,Int}) =
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
Returns true if the commute `u` -> `w` -> `v` obeys `distparam` and `angleparam`.

Uses `gridtype` to determine the distances and angle between `u,w` and `w,v`.
"""
function validtransfer(
        np::TN.TransitNetworkProblem, 
        u::Int,
        v::Int,
        w::Int,
        angleparam::Float64,
        distparam::Float64,
        gridtype::Symbol
    )
    dir1 = dir(np, u, w, gridtype)
    dir2 = dir(np, w, v, gridtype)
    d_uw = edgecost(np, u, w, gridtype)
    d_wv = edgecost(np, w, v, gridtype)
    d_uv = edgecost(np, u, v, gridtype)
    (d_uw + d_wv <= distparam*d_uv) && 
        (dot(dir1,dir2) / norm(dir1) / norm(dir2) >= angleparam)
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
