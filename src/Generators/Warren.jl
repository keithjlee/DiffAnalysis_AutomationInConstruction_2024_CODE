struct Warren2D <: AbstractGenerator
    model::TrussModel
    n::Integer
    dx::Real
    dy::Real
    section::Asap.AbstractSection
    type::Symbol
end

"""
    generatewarren2d(n::Integer,...)

Generate a 2D warren truss in the XY plane.

Required inputs:
- `n::Integer` Number of bays
- `dx::Real` Bay span
- `dy::Real` truss depth
- `section::Asap.AbstractSection` cross section of elements

Default inputs:
- `type::Symbol = :arch` :arch = long chord at bottom; :catenary = long chord at top
"""
function Warren2D(nx::Integer,
        dx::Real,
        dy::Real,
        section::Asap.AbstractSection;
        load = [0., -1., 0.],
        type = :arch)

    @assert type == :arch || type == :catenary "type must be :arch or :catenary"

    #counters
    count = 1
    longids = Vector{Int64}()

    #node collector
    nodes = Vector{TrussNode}()

    #generate longer chord first
    if type == :arch
        longid = :bottomchord
        shortid = :topchord
        y = dy
    else
        longid = :topchord
        shortid = :bottomchord
        y = -dy
    end

    for i = 1:nx+1
        xposition = dx * (i - 1)

        node = TrussNode([xposition, 0., 0.], :free)
        if i == 1
            node.dof = [false, false, false]
            node.id = :pin
        elseif i == nx+1
            node.dof = [true, false, false]
            node.id = :roller
        else
            node.id = longid
        end

        push!(nodes, node)
        push!(longids, count)

        count += 1
    end

    #generate shorter chord
    shortids = Vector{Int64}()
    x0 = dx / 2

    for i = 1:nx
        xposition = x0 + dx * (i - 1)

        node = TrussNode([xposition, y, 0.], :free)
        node.id = shortid

        push!(nodes, node)
        push!(shortids, count)
        count += 1
    end

    #elements
    elements = Vector{TrussElement}()
    
    #long chords
    for i = 1:nx
        element = TrussElement(nodes[longids[i:i+1]]..., section)
        element.id = longid

        push!(elements, element)
    end

    #short chords
    for i = 1:nx-1
        element = TrussElement(nodes[shortids[i:i+1]]..., section)
        element.id = shortid

        push!(elements, element)
    end

    #webs
    for i = 1:nx
        element = TrussElement(nodes[[longids[i], shortids[i]]]..., section)
        element.id = :web
        push!(elements, element)

        element = TrussElement(nodes[[shortids[i], longids[i+1]]]..., section)
        element.id = :web
        push!(elements, element)
    end

    #dummy load
    loads = [NodeForce(n, load) for n in nodes[longid]]

    #assemble and solve
    model = TrussModel(nodes, elements, loads)
    planarize!(model)
    Asap.solve!(model)

    #collect data
    truss = Warren2D(model, nx, dx, dy, section, type)

    #output
    return truss
end

function left_right_mid(ids::Vector{Int64})
    i_left = Int64[]
    i_right = Int64[]
    i_mid = Int64[]

    n_id = length(ids)
    if iseven(n_id)
        n_half = Int(n_id / 2)
        i_left = ids[1:n_half]
        i_right = reverse(ids[n_half+1:end])
    else
        n_mid = Int(ceil(n_id / 2))
        i_left = ids[1:n_mid-1]
        i_right = reverse(ids[n_mid+1:end])
        i_mid = [ids[n_mid]]
    end

    return i_left, i_right, i_mid
end

export find_free_node_pairs
function find_free_node_pairs(gen::Warren2D)

    #bottom chord
    id_bottom = getproperty.(gen.model.nodes[:bottomchord], :nodeID)
    i_bottom_left, i_bottom_right, i_bottom_mid = left_right_mid(id_bottom)

    #top chord
    id_top = getproperty.(gen.model.nodes[:topchord], :nodeID)
    i_top_left, i_top_right, i_top_mid = left_right_mid(id_top)

    # results
    return (
        i_bottom_left = i_bottom_left,
        i_bottom_mid = i_bottom_mid,
        i_bottom_right = i_bottom_right,
        i_top_left = i_top_left,
        i_top_mid = i_top_mid,
        i_top_right = i_top_right
    )

end

export find_element_pairs
function find_element_pairs(gen::Warren2D; tol = 1e-4)

    L = gen.n * gen.dx
    xoffset = L / 2

    midpoints = midpoint.(gen.model.elements)

    xy = [(getindex.(midpoints, 1) .- xoffset) getindex.(midpoints, 2)]

    i = collect(1:gen.model.nElements)

    #there will be only one element that is centered at 0 in x axis
    imid = argmin(abs.(xy[:, 1]))

    i_left = findall(xy[:, 1] .< -tol)
    i_right = findall(xy[:, 1] .> tol)

    @assert length(i_left) == length(i_right)

    xy_left = xy[i_left, :]
    xy_right = xy[i_right, :] .* [-1 1]

    distmat = [norm(e1 - e2) for e1 in eachrow(xy_left), e2 in eachrow(xy_right)]

    ileft = i[i_left]
    iright = i_right[argmin.(eachrow(distmat))]

    return (
        i_left = ileft,
        i_mid = imid,
        i_right = iright
    )
end