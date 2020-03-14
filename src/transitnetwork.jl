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
