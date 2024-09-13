abstract type AbstractGenerator end

export Frame
export SpaceFrame
export Warren2D
export SpaceFrameBeam
export BakerTruss
export TrussFrame
export GridNetwork
export GridFrame
export SimpleTransmissionTower
include("Frame.jl")
include("Spaceframe.jl")
include("Warren.jl")
include("SpaceframeBeam.jl")
include("BakerTruss.jl")
include("TrussFrame.jl")
include("GridNetwork.jl")
include("GridFrame.jl")
include("SimpleTransmissionTower.jl")

export XGroundStructure
export DenseGroundStructure
export BoundedGroundStructure
include("GroundStructure.jl")

export to_truss
export to_frame
include("GroundStructureConversion.jl")