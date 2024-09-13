abstract type AbstractGeo end

include("TrussGeometry.jl")
include("ModelGeometry.jl")
include("NetworkGeometry.jl")

function Geo end

Geo(model::TrussModel) = TrussGeo(model)
Geo(model::Model) = ModelGeo(model)
Geo(network::Network) = NetworkGeo(network)