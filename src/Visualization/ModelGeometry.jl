struct ModelGeo <: AbstractGeo
    nodes::Vector{Vector{Float64}}
    nodes_xy::Vector{Vector{Float64}}
    disp::Vector{Vector{Float64}}
    disp_xy::Vector{Vector{Float64}}
    indices::Vector{Vector{Int64}}
    indices_flat::Vector{Int64}
    P::Vector{Float64}
    max_abs_P::Float64
    Vy::Vector{Float64}
    max_abs_Vy::Float64
    Vz::Vector{Float64}
    max_abs_Vz::Float64
    Tx::Vector{Float64}
    max_abs_Tx::Float64
    My::Vector{Float64}
    max_abs_My::Float64
    Mz::Vector{Float64}
    max_abs_Mz::Float64
    areas::Vector{Float64}
    max_area::Float64
    Ix::Vector{Float64}
    max_Ix::Float64
    Iy::Vector{Float64}
    max_Iy::Float64
    J::Vector{Float64}
    max_J::Float64
    lengths::Vector{Float64}
    element_vectors::Vector{Vector{Float64}}
    element_vectors_xy::Vector{Vector{Float64}}

    function ModelGeo(model::Model)

        nodes = getproperty.(model.nodes, :position)
        nodes_xy = [node[1:2] for node in nodes]

        disp = [node.displacement[1:3] for node in model.nodes]
        disp_xy = [d[1:2] for d in disp]

        indices = Asap.nodeids.(model.elements)
        indices_flat = vcat(indices...)

        element_forces = getproperty.(model.elements, :forces)

        P = vcat([forces[[1,7]] .* [-1, 1] for forces in element_forces]...)
        max_abs_P = maximum(abs.(P))

        Vy = vcat([forces[[2,8]] .* [-1, 1] for forces in element_forces]...)
        max_abs_Vy = maximum(abs.(Vy))

        Vz = vcat([forces[[3,9]] .* [-1, 1] for forces in element_forces]...)
        max_abs_Vz = maximum(abs.(Vz))

        Tx = vcat([forces[[4,10]] .* [-1, 1] for forces in element_forces]...)
        max_abs_Tx = maximum(abs.(Tx))

        My = vcat([forces[[5,11]] .* [-1, 1] for forces in element_forces]...)
        max_abs_My = maximum(abs.(My))

        Mz = vcat([forces[[6,12]] .* [-1, 1] for forces in element_forces]...)
        max_abs_Mz = maximum(abs.(Mz))

        sections = getproperty.(model.elements, :section)

        areas = getproperty.(sections, :A)
        max_area = maximum(areas)

        Ix = getproperty.(sections, :Ix)
        max_Ix = maximum(Ix)

        Iy = getproperty.(sections, :Iy)
        max_Iy = maximum(Iy)

        J = getproperty.(sections, :J)
        max_J = maximum(J)

        element_vectors = Asap.local_x.(model.elements)
        element_vectors_xy = [evec[1:2] for evec in element_vectors]

        return new(
            nodes,
            nodes_xy,
            disp,
            disp_xy,
            indices,
            indices_flat,
            P,
            max_abs_P,
            Vy,
            max_abs_Vy,
            Vz,
            max_abs_Vz,
            Tx,
            max_abs_Tx,
            My,
            max_abs_My,
            Mz,
            max_abs_Mz,
            areas,
            max_area,
            Ix,
            max_Ix,
            Iy,
            max_Iy,
            J,
            max_J,
            getproperty.(model.elements, :length),
            element_vectors,
            element_vectors_xy
        )
    end

end