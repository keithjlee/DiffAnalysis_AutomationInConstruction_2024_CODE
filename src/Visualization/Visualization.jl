include("Geometry.jl")
export Geo, ModelGeo, TrussGeo, NetworkGeo

include("Colours.jl")
export kjl_pink, kjl_blue, kjl_green, kjl_orange, kjl_gray, kjl_tan, kjl_turquoise
export caitlin_blue, caitlin_pink
export pink2blue
export white2black
export white2blue
export white2pink

include("Themes.jl")
export aic

halfwidth(ratio::Real = 1.0; factor = 2) = Int.(round.(factor .* 255 .* (1, ratio)))
export halfwidth

fullwidth(ratio::Real = 1.0; factor = 2) = Int.(round.(factor .* 539 .* (1, ratio)))
export fullwidth

include("Utils.jl")
export style1!
export discretize
export labelscale!
export gridtoggle!
export simplifyspines!
export mirrorticks!
export alignticks!
export tickstoggle!