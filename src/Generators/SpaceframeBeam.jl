struct SpaceFrameBeam <: AbstractGenerator
    model::TrussModel
    nx::Integer
    dx::Real
    dy::Real
    dz::Real
    section::Asap.AbstractSection
    
    function SpaceFrameBeam(n::Integer, dx::Real, dy::Real, dz::Real, section::Asap.AbstractSection, load = [0., 0., -10.])

        @assert n % 2 == 0 "n must be even"

        #x positions
        x_positions = collect(0:dx:n*dx)

        #y positions
        y_bottom = zero(x_positions)
        y_top1 =  fill(dy/2, length(x_positions))
        y_top2 = fill(-dy/2, length(x_positions))

        #z positions
        z_bottom = zero(x_positions)
        z_top = fill(dz, length(x_positions))

        #make bottom nodes
        bottom_nodes = [TrussNode([x, y, z], :free, :bottom) for (x,y,z) in zip(x_positions, y_bottom, z_bottom)]

        #top nodes 1
        top1_nodes = [TrussNode([x,y,z], :free, :top1) for (x,y,z) in zip(x_positions, y_top1, z_top)]
        first(top1_nodes).id = :pin
        fixnode!(first(top1_nodes), :pinned)
        last(top1_nodes).id = :roller
        fixnode!(last(top1_nodes), :xfree)

        #top nodes 2
        top2_nodes = [TrussNode([x,y,z], :free, :top2) for (x,y,z) in zip(x_positions, y_top2, z_top)]
        first(top2_nodes).id = :pin
        fixnode!(first(top2_nodes), :pinned)
        last(top2_nodes).id = :roller
        fixnode!(last(top2_nodes), :xfree)

        #make bottom elements
        bottom_elements = [TrussElement(bottom_nodes[i], bottom_nodes[i+1], section, :bottom) for i = 2:length(bottom_nodes)-2]

        #make top elements
        top1_elements = [TrussElement(top1_nodes[i], top1_nodes[i+1], section, :top1) for i = 1:n]
        top2_elements = [TrussElement(top2_nodes[i], top2_nodes[i+1], section, :top2) for i = 1:n]

        #top struts
        strut_elements = [TrussElement(node1, node2, section, :strut) for (node1, node2) in zip(top1_nodes, top2_nodes)]

        #web elements
        web_elements = Vector{TrussElement}()

        for i = 1:Int(n/2)

            push!(web_elements, TrussElement(top1_nodes[i], bottom_nodes[i+1], section, :web))
            push!(web_elements, TrussElement(top2_nodes[i], bottom_nodes[i+1], section, :web))

        end

        for i = Int(n/2)+1:n

            push!(web_elements, TrussElement(bottom_nodes[i], top1_nodes[i+1], section, :web))
            push!(web_elements, TrussElement(bottom_nodes[i], top2_nodes[i+1], section, :web))

        end

        for i = 2:n
            push!(web_elements, TrussElement(top1_nodes[i], bottom_nodes[i], section, :web))
            push!(web_elements, TrussElement(top2_nodes[i], bottom_nodes[i], section, :web))
        end

        nodes = [bottom_nodes[2:end-1]; top1_nodes; top2_nodes]
        elements = [bottom_elements; top1_elements; top2_elements; strut_elements; web_elements]
        loads = [NodeForce(node, load) for node in nodes[:bottom]]

        model = TrussModel(nodes, elements, loads)
        Asap.solve!(model)

        return new(
            model,
            n,
            dx,
            dy,
            dz,
            section
        )
    end
end

