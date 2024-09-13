struct HalfDome{M,S} <: AbstractGenerator
    model::M
    radius::Float64
    n_circ::Integer
    n_rings::Integer
    section::S
end

function HalfDome(radius::Float64, n_circ::Integer, n_rings::Integer, section::Section; load = [0., 0., -10.], support = :pinned)

    horizontal_angle_range = range(0, 2pi, n_circ+1)[1:end-1]
    vertical_angle_range = range(0, pi/2, n_rings+1)[1:end-1]

    horizontal_offset_angle = horizontal_angle_range[2]

    #make nodes
    node_collector = Vector{Vector{Node}}()

    for i = 1:n_rings

        ϕ = vertical_angle_range[i]

        #checks
        offset = i % 2 == 0 ? horizontal_offset_angle : 0.
        id = i % 2 == 0 ? :free : base
        dof = i % 2 == 0 ? :free : support

        xy_factor = cos(ϕ) * radius
        z = sin(ϕ) * radius

        nodes = [
            Node([xy_factor * cos(θ + offset), xy_factor * sin(θ + offset), z], dof, id) for θ in horizontal_angle_range
        ]

        push!(node_collector, nodes)
    end

    #last node
    push!(node_collector, [Node(0., 0., radius), :free, :top])

    
end