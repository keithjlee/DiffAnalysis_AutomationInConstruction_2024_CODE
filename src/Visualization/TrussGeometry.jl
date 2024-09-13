"""
    TrussGeo

Utility structure for visualizing truss models
"""
struct TrussGeo <: AbstractGeo
    nodes::Vector{Vector{Float64}}
    nodes_xy::Vector{Vector{Float64}}
    disp::Vector{Vector{Float64}}
    disp_xy::Vector{Vector{Float64}}
    indices::Vector{Vector{Int64}}
    indices_flat::Vector{Int64}
    forces::Vector{Float64}
    max_abs_force::Float64
    areas::Vector{Float64}
    max_area::Float64
    lengths::Vector{Float64}
    element_vectors::Vector{Vector{Float64}}
    element_vectors_xy::Vector{Vector{Float64}}
    load_positions::Vector{Vector{Float64}}
    load_positions_xy::Vector{Vector{Float64}}
    load_vectors::Vector{Vector{Float64}}
    load_vectors_xy::Vector{Vector{Float64}}

    function TrussGeo(model::TrussModel)

        nodes = getproperty.(model.nodes, :position)
        nodes_xy = [node[1:2] for node in nodes]

        disp = getproperty.(model.nodes, :displacement)
        disp_xy = [d[1:2] for d in disp]

        indices = Asap.nodeids.(model.elements)
        indices_flat = vcat(indices...)

        forces = getindex.(getproperty.(model.elements, :forces), 2)
        max_abs_force = maximum(abs.(forces))

        areas = getproperty.(getproperty.(model.elements, :section), :A)
        max_area = maximum(areas)

        load_positions = getproperty.(getproperty.(model.loads, :node), :position)
        load_positions_xy = [load[1:2] for load in load_positions]

        load_vectors = getproperty.(model.loads, :value)
        load_vectors_xy = [load[1:2] for load in load_vectors]

        element_vectors = Asap.local_x.(model.elements)
        element_vectors_xy = [evec[1:2] for evec in element_vectors]


        return new(
            nodes,
            nodes_xy,
            disp,
            disp_xy,
            indices,
            indices_flat,
            forces,
            max_abs_force,
            areas,
            max_area,
            getproperty.(model.elements, :length),
            element_vectors,
            element_vectors_xy,
            load_positions,
            load_positions_xy,
            load_vectors,
            load_vectors_xy
        )

    end

end