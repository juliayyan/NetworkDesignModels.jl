"""
TransitNetwork

* `nstations`: number of stations in network
* `lines`: lines in the network (which can be used to start a master problem)
* `xfrstns`: for each commute (u, v), a list of transfer stations w that could be used to 
*            transfer u -> w -> v.
* `odmatrix`: demand for commute (u, v)
* `latlon`: (latlon[u, 1], latlon[u, 2]) provide the coordinates for station u
* `gridtype`: `:euclidan` if `latlon` should be interpreted as Euclidean coordinates and `:latlong" otherwise
* `dists`: (if not specified, then computed using haversine or euclidean distance)
*           dists[u,v] gives the edge distance for edge (u,v)
* `spdists`: (if not specified, then computed using haversine or euclidean distance)
*            spdists[u,v] gives the shortest path distance from u to v
"""
mutable struct TransitNetwork
    nstations::Int
    lines::Vector{Vector{Int}}
    xfrstns::Dict{Tuple{Int,Int},Vector{Int}}
    odmatrix::Array{Float64,2}
    latlon::Array{Float64,2}
    gridtype::Symbol
    dists
    spdists
end

function Base.show(io::IO, np::TransitNetwork; offset::String="")
    println(io, offset, string(typeof(np)))
    println(io, offset, "    * nstations: $(np.nstations)")
end
