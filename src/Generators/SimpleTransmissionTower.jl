
function SimpleTransmissionTower(dx::Float64, dy::Float64, ny::Int64, doffset::Float64, sec::TrussSection; Pcable = [0., -1., 0.])


    x_col1 = -dx / 2
    x_col2 = dx / 2

    yvals = range(0, dy*ny, ny+1)
    yvals2 = range(dy/2, dy*ny-dy/2, ny)

    #left nodes
    nodes_left = [TrussNode([x_col1, y, 0.], :free, :col1) for y in yvals]

    #right nodes
    nodes_right = [TrussNode([x_col2, y, 0.], :free, :col2) for y in yvals]

    #middle nodes
    nodes_middle = [TrussNode([0., y, 0.], :free, :middle) for y in yvals2]

    #offset nodes
    node_left_offset = TrussNode([x_col1 - doffset, yvals[end-1], 0.], :free, :leftoffset)
    node_right_offset = TrussNode([x_col2 + doffset, yvals[end-1], 0.], :free, :rightoffset)

    #offset midpoints
    node_left_offset_midpoint = TrussNode([x_col1 - 0.5doffset, yvals[end-1] - dy/2, 0.], :free, :leftoffsetmid)
    node_right_offset_midpoint = TrussNode([x_col2 + 0.5doffset, yvals[end-1] - dy/2, 0.], :free, :rightoffsetmid)

    #elements

    #left vertical
    elements_left = [TrussElement(nodes_left[i], nodes_left[i+1], sec, :col1) for i = 1:ny]
    
    #right vertical
    elements_right = [TrussElement(nodes_right[i], nodes_right[i+1], sec, :col2) for i = 1:ny]

    #horizontal
    elements_horizontal = [TrussElement(nl, nr, sec, :horizontal) for (nl, nr) in zip(nodes_left[2:end], nodes_right[2:end])]

    #web left 1
    elements_web_left_1 = [TrussElement(nl, nm, sec, :webleft) for (nl, nm) in zip(nodes_left[1:end-1], nodes_middle)]
    
    #web left 2
    elements_web_left_2 = [TrussElement(nl, nm, sec, :webleft) for (nl, nm) in zip(nodes_left[2:end], nodes_middle)]

    #web right 1
    elements_web_right_1 = [TrussElement(nr, nm, sec, :webright) for (nr, nm) in zip(nodes_right[1:end-1], nodes_middle)]
    
    #web right 2
    elements_web_right_2 = [TrussElement(nr, nm, sec, :webright) for (nr, nm) in zip(nodes_right[2:end], nodes_middle)]

    #left offset
    elements_left_offset = [
        TrussElement(node_left_offset, nodes_left[end-1], sec, :leftoffset),
        TrussElement(node_left_offset, node_left_offset_midpoint, sec, :leftoffset),
        TrussElement(node_left_offset_midpoint, nodes_left[end-2], sec, :leftoffset),
        TrussElement(node_left_offset_midpoint, nodes_left[end-1], sec, :leftoffset)
    ]

    #right offset
    elements_right_offset = [
        TrussElement(node_right_offset, nodes_right[end-1], sec, :rightoffset),
        TrussElement(node_right_offset, node_right_offset_midpoint, sec, :rightoffset),
        TrussElement(node_right_offset_midpoint, nodes_right[end-2], sec, :rightoffset),
        TrussElement(node_right_offset_midpoint, nodes_right[end-1], sec, :rightoffset)
    ]

    #supports
    nodes_left[1].id = :leftsupport
    nodes_left[1].dof = [false, false, false]
    nodes_right[1].id = :rightsupport
    nodes_right[1].dof = [false, false, false]

    #collect
    nodes = [nodes_left; nodes_right; nodes_middle; node_left_offset; node_left_offset_midpoint; node_right_offset; node_right_offset_midpoint]
    elements = [elements_left; elements_right; elements_horizontal; elements_web_left_1; elements_web_left_2; elements_web_right_1; elements_web_right_2; elements_left_offset; elements_right_offset]

    #loads
    loads = [NodeForce(node_left_offset, Pcable), NodeForce(node_right_offset, Pcable)]

    #assemble
    model = TrussModel(nodes, elements, loads)

    planarize!(model)
    Asap.solve!(model)
    

    return model
end