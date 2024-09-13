"""
    Frame <: AbstractGenerator

A multistory frame structure with joist-beam-column hierarchy.

# Fields
- `model` Asap model
- `nx` number of bays in x direction
- `dx` bay spacing in x direction
- `ny` number of bays in y direction
- `dy` bay spacing in y direction
- `dz` column height
- `joistSpacing` spacing of joists
- `columnSection` Asap section used for columns
- `primarySection` Asap section used for primary spans
- `joistSection` Asap section used for joists
- `braceSection` Asap section used for braces
- `columnRelease` End release for columns
- `primaryRelease` End release used for primary spans
- `joistRelease` End release used for joists
- `braceRelease` End release used for braces
- `columnPsi` LCS rotation for columns
- `primaryPsi` LCS rotation for primary spans
- `joistPsi` LCS rotation for joists
- `Base` Base point of Frame generator
- `iExteriorXnodes` Indices of nodes that are parallel to the global X axis, are on the exterior ends of the frame, and coincide with primary-column intersections
- `iExteriorYnodes` Indices of nodes that are parallel to the global Y axis, are on the exterior ends of the frame, and coincide with joist-column intersections
- `iExteriorXjoists` Indices of elements that are parallel to the global Y axis and are on the exterior ends of the frame
- `iExteriorYprimaries` indices of elements that are parallel to the global X axis and are on the exterior ends of the frame
"""
struct Frame <: AbstractGenerator
    model::Model
    nx::Integer
    dx::Real
    ny::Integer
    dy::Real
    nz::Integer
    dz::Real
    joistSpacing::Real
    columnSection::Section
    primarySection::Section
    joistSection::Section
    braceSection::Section
    columnRelease::Symbol
    primaryRelease::Symbol
    joistRelease::Symbol
    braceRelease::Symbol
    columnPsi::Real
    primaryPsi::Real
    joistPsi::Real
    base::Vector{<:Real}
    iExteriorXnodes::Vector{Integer}
    iExteriorYnodes::Vector{Integer}
    iExteriorXjoists::Vector{Integer}
    iExteriorYprimaries::Vector{Integer}
end

"""
    Frame(nx, dx, ny, dy, nz, dz, joistSpacing, columnSection, primarySection, joistSection, braceSection)

Generate a 3D frame model.

# Arguments
- `nx::Integer` Number of bays in primary span direction
- `dx::Real` Primary span length
- `ny::Integer` Number of bays in secondary span direction
- `dy::Real` Secondary span length
- `nz::Integer` Number of stories
- `dz::Real` Story height
- `joistSpacing::Real` Approx. center spacing between secondary spans
- `columnSection::Section` Section of column elements
- `primarySection::Section` Section of primary beams
- `joistSection::Section` Section of secondary beams
- `braceSection::Section` Section for lateral bracing

# Optional inputs
- `diaphragm::Bool = false` Include floor diaphragm (as X braces between columns)
- `diaphragmSection::Section = nothing` Element section for diaphragm members (defaults to braceSection)
- `columnRelease::Symbol = :fixedfixed` DOF end release for column elements
- `primaryRelease::Symbol = :fixedfixed` DOF end release for primary elements
- `joistRelease::Symbol = :joist` DOF end release for joist elements
- `braceRelease::Symbol = :freefree` DOF end release for braces
- `columnPsi::Real = 0` Angle of roll Ψ for column LCS
- `primaryPsi::Real = π/2` Angle of roll Ψ for primary beam LCS
- `joistPsi::Real = π/2` Angle of roll Ψ for secondary beam LCS
- `base::Vector{Real} = [0, 0, 0]` base point for frame grid generation 
"""
function Frame(nx::Integer,
        dx::Real,
        ny::Integer,
        dy::Real,
        nz::Integer,
        dz::Real,
        joistSpacing::Real,
        columnSection::Section,
        primarySection::Section,
        joistSection::Section,
        braceSection::Section;
        diaphragm = false,
        diaphragmSection = nothing,
        columnRelease = :fixedfixed,
        primaryRelease = :fixedfixed,
        joistRelease = :joist,
        braceRelease = :freefree,
        columnPsi = 0.,
        primaryPsi = pi/2,
        joistPsi = pi/2,
        base = [0., 0., 0.]
        )


    ########
    # nodes
    ########

    xoffset = [dx, 0., 0.]
    yoffset = [0., dy, 0.]
    zoffset = [0., 0., dz]

    #generate
    nodes = [Node(base + xoffset * i + yoffset * j + zoffset * k, :pinned) for i in 0:nx, j in 0:ny, k in 0:nz]

    #release non-ground nodes
    for node in nodes
        if last(node.position) > last(base)
            fixnode!(node, :free)
        end
    end

    ######
    # columns
    #######
    #make columns
    columns = Vector{Element}()
    for i in 1:nx+1
        for j in 1:ny+1
            for k in 1:nz
                el = Element(nodes[i,j,k], nodes[i,j,k+1], columnSection, columnRelease)
                # el.Ψ = 0
                el.Ψ = columnPsi
                el.id = :column
                push!(columns, el)
            end
        end
    end


    #####
    # primary beams
    #####
    primaries = Vector{Element}()

    for k in 2:nz+1
        for i in 1:nx
            for j in 1:ny+1
                el = Element(nodes[i,j,k], nodes[i+1,j,k], primarySection, primaryRelease)
                el.id = :primary
                el.Ψ = primaryPsi
                push!(primaries, el)
            end
        end
    end

    ######
    # joists
    ######

    primreshaped = reshape(primaries, ny+1, nx, nz)

    secondaries = Vector{Union{BridgeElement, Element}}()
    njoists = Int(round(dx / joistSpacing))

    #main bridge elements
    for k in 1:nz
        for j = 1:nx
            for i = 1:ny
                bridges = [BridgeElement(primreshaped[i,j,k], x, primreshaped[i+1,j,k], x, joistSection, joistRelease) for x in range(0,1,njoists)[2:end-1]]

                for bridge in bridges
                    bridge.id = :joist
                    bridge.Ψ = joistPsi
                end

                push!(secondaries, bridges...)

            end
        end
    end

    #between column nodes
    for i in 1:nx+1
        for j in 1:ny
            for k in 2:nz+1
                bridge = Element(nodes[i,j,k], nodes[i,j+1,k], joistSection, joistRelease)
                # bridge.Ψ = 0.
                bridge.id = :joist
                bridge.Ψ = joistPsi
                push!(secondaries, bridge)
            end
        end
    end

    #########
    #braces
    #########

    braces = Vector{Element}()
    i = 1
    for j = [1, ny+1]
        for k = 1:nz
            brace1 = Element(nodes[i,j,k], nodes[i+1,j,k+1], braceSection, braceRelease)
            brace2 = Element(nodes[i+1,j,k], nodes[i,j,k+1], braceSection, braceRelease)

            brace1.id = brace2.id = :brace
            push!(braces, brace1, brace2)
        end
    end

    j = 1
    for i = [1, nx+1]
        for k = 1:nz
            brace1 = Element(nodes[i,j,k], nodes[i, j+1, k+1], braceSection, braceRelease)
            brace2 = Element(nodes[i,j+1,k], nodes[i,j,k+1], braceSection, braceRelease)

            brace1.id = brace2.id = :brace
            push!(braces, brace1, brace2)
        end
    end

    ######
    #diaphragm
    if diaphragm

        dsection = isnothing(diaphragmSection) ? braceSection : diaphragmSection

        for k = 2:nz+1
            for i = 1:nx
                for j = 1:ny
                    diaph1 = Element(nodes[i, j, k], nodes[i+1, j+1, k], dsection, :freefree)
                    diaph2 = Element(nodes[i+1, j, k], nodes[i, j+1, k], dsection, :freefree)

                    diaph1.id = diaph2.id = :diaphragm
                    push!(braces, diaph1, diaph2)
                end
            end
        end
    end


    ######
    # dummy load and assembly
    ######

    loads = [LineLoad(j, [0., 0., -1.]) for j in secondaries]
    flatnodes = vec(nodes)
    elements = [columns; primaries; secondaries; braces]

    model = Model(flatnodes, elements, loads)
    Asap.solve!(model)


    #extract node/element indices
    iExteriorXnodes = [getproperty.(vec(nodes[:,1,:]), :nodeID); getproperty.(vec(nodes[:,end,:]), :nodeID)]
    iExteriorYnodes = [getproperty.(vec(nodes[1,:,:]), :nodeID); getproperty.(vec(nodes[end,:,:]), :nodeID)]

    iExteriorXjoists = Vector{Int64}()
    iExteriorYprimary = Vector{Int64}()

    xExtrema = [first(base), first(base) + dx * nx]
    yExtrema = [base[2], base[2] + dy * ny]
    for (i,element) in enumerate(model.elements)
        if element.id == :joist
            x = element.nodeStart.position[1]
            if minimum(abs.(xExtrema .- x)) <= model.tol
                push!(iExteriorXjoists, i)
            end
        elseif element.id == :primary
            y = element.nodeStart.position[2]
            if minimum(abs.(yExtrema .- y)) <= model.tol
                push!(iExteriorYprimary, i)
            end
        end
    end

    frame = Frame(model,
        nx,
        dx,
        ny,
        dy,
        nz,
        dz,
        joistSpacing,
        columnSection,
        primarySection,
        joistSection,
        braceSection,
        columnRelease,
        primaryRelease,
        joistRelease,
        braceRelease,
        columnPsi,
        primaryPsi,
        joistPsi,
        base,
        iExteriorXnodes,
        iExteriorYnodes,
        iExteriorXjoists,
        iExteriorYprimary)

    return frame;
end