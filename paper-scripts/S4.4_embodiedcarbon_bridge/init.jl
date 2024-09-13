using DiffAnalysis_AIC2024
set_theme!(aic)
using Nonconvex, NonconvexNLopt

datadir = joinpath(@__DIR__, "data/"); isdir(datadir) || mkdir(datadir)
figdir = joinpath(@__DIR__, "figures/"); isdir(figdir) || mkdir(figdir)

# material properties and initial sections
begin
    E_glulam = 13.1e6 #kN/m²
    ρ_glulam = 560. #kg/m³ [Nordic Lam+]
    ECC_glulam = ρ_glulam * 0.512 #kg CO₂e / kg [0.512 kg/kg taken from ICE DB (glulam, no offset)]
    b0_glulam = .4
    d0_glulam = 1.25
    A0_glulam = a_box(b0_glulam, d0_glulam)
    ft_glulam = 20.4e3 # kN/m²
    fc_glulam = 33e3 # kN/m²

    init_glulam_section = TrussSection(A0_glulam, E_glulam)

    E_steel = 200e6 #kN/m²
    ρ_steel = 7800. #kg/m³
    ECC_steel = ρ_steel * 1.55 #kg CO₂e / kg [1.55 kg/kg taken from ICE DB (steel section)]
    b0_steel = .2
    d0_steel = .7
    t0_steel = .02
    A0_steel = a_box(b0_steel, d0_steel) - a_box(b0_steel - 2t0_steel, d0_steel - 2t0_steel)
    fy_steel = 350e3 #MPa

    init_steel_section = TrussSection(A0_steel, E_steel)
end

# design parameters

L = 120. # total span
L_sidespine = 30. # side span (approx)
support_vertical_offset = 12. # height of bridge above ground

nx = 39 # number of bays
dy = 5. # depth 
P = 150. # deck load

dx = L / nx

dispmax = (L - 2L_sidespine) / 400

# make structure
begin
    # bottom chord nodes
    xvals = [dx * (i-1) for i = 1:nx+1]
    bottom_nodes = [TrussNode([x, 0., 0.], :free, :bottomchord) for x in xvals]

    first(bottom_nodes).dof = [true, false, false]
    last(bottom_nodes).dof = [true, false ,false]

    # supports
    i_support_offset = Int(div(L_sidespine, dx))

    i_support_left = i_support_offset:i_support_offset+1
    x_support_left = xvals[i_support_left]
    i_support_right = nx+1-i_support_offset:nx+2-i_support_offset
    x_support_right = xvals[i_support_right]

    #left supports
    left_support_nodes = [TrussNode([x, -support_vertical_offset, 0.], :pinned, :leftsupport) for x in x_support_left]
    left_support_mid = TrussNode([mean(x_support_left), -.5support_vertical_offset, 0.], :free, :leftsupportmid)

    #right supports
    right_support_nodes = [TrussNode([x, -support_vertical_offset, 0.], :pinned, :rightsupport) for x in x_support_right]
    right_support_mid = TrussNode([mean(x_support_right), -.5support_vertical_offset, 0.], :free, :rightsupportmid)

    # middle nodes
    middle_nodes = [
        TrussNode([-dx, dy/2, 0.], :xfree, :leftspanend);
        [TrussNode([x + dx/2, dy/2, 0.], :free, :middle) for x in xvals[1:end-1]];
        TrussNode([L + dx, dy/2, 0.], :xfree, :rightspanend)
    ]

    # top nodes
    top_nodes = [TrussNode([x, dy, 0.], :free, :topchord) for x in xvals]

    # chord
    bottom_elements = [TrussElement(bottom_nodes[i], bottom_nodes[i+1], init_glulam_section, :bottomchord) for i = 1:nx]

    top_elements  = [TrussElement(top_nodes[i], top_nodes[i+1], init_glulam_section, :topchord) for i = 1:nx]

    # webs
    bottom_diagonals = [
        [TrussElement(n1, n2, init_glulam_section, :bottomweb) for (n1, n2) in zip(bottom_nodes[1:end-1], middle_nodes[2:end-1])];
        [TrussElement(n1, n2, init_glulam_section, :bottomweb) for (n1, n2) in zip(bottom_nodes[2:end], middle_nodes[2:end-1])]
    ]

    top_diagonals = [
        [TrussElement(n1, n2, init_glulam_section, :topweb) for (n1, n2) in zip(top_nodes[1:end-1], middle_nodes[2:end-1])];
        [TrussElement(n1, n2, init_glulam_section, :topweb) for (n1, n2) in zip(top_nodes[2:end], middle_nodes[2:end-1])]
    ]

    # verticals
    verticals = [TrussElement(n1, n2, init_glulam_section, :strut) for (n1, n2) in zip(bottom_nodes, top_nodes)]

    # support elements
    support_elements = [
        [TrussElement(n1, n2, init_glulam_section, :support) for (n1, n2) in zip(left_support_nodes, bottom_nodes[i_support_left])];
        [TrussElement(n1, n2, init_glulam_section, :support) for (n1, n2) in zip(right_support_nodes, bottom_nodes[i_support_right])];
        TrussElement(first(left_support_nodes), bottom_nodes[first(i_support_left)-1], init_glulam_section, :support);
        TrussElement(last(left_support_nodes), bottom_nodes[last(i_support_left)+1], init_glulam_section, :support);
        TrussElement(first(right_support_nodes), bottom_nodes[first(i_support_right)-1], init_glulam_section, :support);
        TrussElement(last(right_support_nodes), bottom_nodes[last(i_support_right)+1], init_glulam_section, :support);
        [TrussElement(node, left_support_mid, init_glulam_section, :support) for node in bottom_nodes[i_support_left]];
        [TrussElement(node, left_support_mid, init_glulam_section, :support) for node in left_support_nodes];
        [TrussElement(node, right_support_mid, init_glulam_section, :support) for node in bottom_nodes[i_support_right]];
        [TrussElement(node, right_support_mid, init_glulam_section, :support) for node in right_support_nodes];
    ]

    #span elements
    deck = [TrussElement(middle_nodes[i], middle_nodes[i+1], init_steel_section, :deck) for i = 1:nx+1]

    nodes = [
        bottom_nodes;
        middle_nodes;
        top_nodes;
        left_support_nodes;
        left_support_mid;
        right_support_nodes;
        right_support_mid
    ]

    loads = [NodeForce(node, [0., -P, 0.]) for node in middle_nodes[2:end-1]]

    elements = [
        bottom_elements;
        bottom_diagonals;
        top_diagonals;
        top_elements;
        verticals;
        support_elements;
        deck
    ]

    model = TrussModel(nodes, elements, loads)
    for node in model.nodes
        node.position .+= [0., support_vertical_offset, 0.]
    end

    planarize!(model)
    Asap.solve!(model)
end

# viz
begin
    geo = Geo(model)
    dfac = Observable(0.)
    pts = @lift(Point2.(geo.nodes_xy .+ $dfac .* geo.disp_xy))
    els = @lift($pts[geo.indices_flat])

    fig = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect())
    hidespines!(ax)
    ylims!(ax, 0., 2dy + support_vertical_offset)

    linesegments!(els)

    sl = Slider(fig[2,1], range = range(0, 100, 100), startvalue = 0)
    on(sl.value) do val
        dfac[] = val
    end

    fig
end

# node pairs for spatial variables
n_top = length(model.nodes[:topchord])
topdiv = div(n_top, 2)
i_top = findall(model.nodes, :topchord)
i_left_top = i_top[1:topdiv]
i_right_top = reverse(i_top[topdiv+1:end])


n_bottom = length(model.nodes[:bottomchord])
bottomdiv = div(n_bottom, 2)
i_bottom = findall(model.nodes, :bottomchord)

i_left_bottom = i_bottom[1:bottomdiv]
i_right_bottom = reverse(i_bottom[bottomdiv+1:end])

i_support_mid_left = left_support_mid.nodeID
i_support_mid_right = right_support_mid.nodeID

i_structure = [i for i in eachindex(model.elements) if model.elements[i].id != :deck]
i_deck = [i for i in eachindex(model.elements) if model.elements[i].id == :deck]

# element indices
i_bottom_elements = findall(model.elements, :bottomchord)
i_top_elements = findall(model.elements, :topchord)
i_topweb_elements = findall(model.elements, :topweb)
i_bottomweb_elements = findall(model.elements, :bottomweb)
i_strut_elements = findall(model.elements, :strut)
i_support_elements = findall(model.elements, :support)

i_deck_nodes = findall(model.nodes, :middle)

material_dict = Dict(:wood => "W", :steel => "S")


#visualization data
ls_indices = vcat(geo.indices[i_structure]...)
ls_deck = vcat(geo.indices[i_deck]...)

static_deck_points = Point2.(geo.nodes_xy[ls_deck])

display(fig)