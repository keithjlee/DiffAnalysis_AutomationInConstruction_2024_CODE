include("Utilities.jl")
export replace_values
export add_values

include("Geometry.jl")

include("Rtruss.jl")
include("Rframe.jl")

include("Ktruss.jl")
include("Kframe.jl")

include("K.jl")

include("Solve.jl")

include("Objective.jl")
export solve_truss
export solve_truss_noadjoint
export compliance

export solve_network
export target

export solve_frame
export solve_frame_noadjoint

include("PostProcessing.jl")
export element_forces
export axial_force
export axial_stress
export updatemodel
export updatenetwork
export element_forces

include("Constraints.jl")
export disp_stress_cstr

include("SectionMaker.jl")
export a_tube
export I_tube
export J_tube
export S_tube

export a_box
export Ix_box
export Iy_box
export J_box
export Sx_box
export Sy_box