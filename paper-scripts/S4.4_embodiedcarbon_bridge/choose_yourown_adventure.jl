date = "09_11_2024"

# choose your material types here
begin
    bottomchord_material = :wood
    topchord_material = :wood
    bottomweb_material = :wood
    topweb_material = :wood
    strut_material = :wood
    support_material = :wood
end

begin
    include("init.jl")
    # collect
    element_material_groups = [bottomchord_material, topchord_material, bottomweb_material, topweb_material, strut_material, support_material]
    element_id_groups = [:bottomchord, :topchord, :bottomweb, :topweb, :strut, :support]

    # update sections (everything is by default wood, so only update if steel is selected)
    for (mat, group) in zip(element_material_groups, element_id_groups)
        if mat == :steel
            for element in model.elements[group]
                element.section = init_steel_section
            end
        end
    end

    #re-solve
    Asap.solve!(model; reprocess = true)
end

#=
VARIABLE DEFINITION
=#

# SPATIAL BOUNDS
begin
    fac = 0.9
    dxmin = -dx / 2 * fac
    dxmax = dx / 2 * fac

    dymin_bottom = -fac * support_vertical_offset
    dymax_bottom = dy / 3

    dymin_top = - dy / 3
    dymax_top = 5.

    dymin_support = -support_vertical_offset/2 * fac
    dymax_support = support_vertical_offset/2 * fac
end

# AREA BOUNDS
begin
    bmin = .2
    bmax = .6
    dmin = .3
    dmax = 3.

    Amin_wood = a_box(bmin, dmin)
    Amax_wood = a_box(bmax, dmax)

    Amin_steel = 1e-3
    Amax_steel = 0.2
end

#apply bounds
begin
    Amin_bottom = bottomchord_material == :steel ? Amin_steel : Amin_wood
    Amax_bottom = bottomchord_material == :steel ? Amax_steel : Amax_wood


    Amin_top = topchord_material == :steel ? Amin_steel : Amin_wood
    Amax_top = topchord_material == :steel ? Amax_steel : Amax_wood

    Amin_bottomweb = bottomweb_material == :steel ? Amin_steel : Amin_wood
    Amax_bottomweb = bottomweb_material == :steel ? Amax_steel : Amax_wood

    Amin_topweb = topweb_material == :steel ? Amin_steel : Amin_wood
    Amax_topweb = topweb_material == :steel ? Amax_steel : Amax_wood

    Amin_strut = strut_material == :steel ? Amin_steel : Amin_wood
    Amax_strut = strut_material == :steel ? Amax_steel : Amax_wood

    Amin_support = support_material == :steel ? Amin_steel : Amin_wood
    Amax_support = support_material == :steel ? Amax_steel : Amax_wood
end

# define variables
begin
# SPATIAL VARIABLES
    x_bottom_parents = [SpatialVariable(node, 0., dxmin, dxmax, :X) for node in model.nodes[i_left_bottom[2:end]]]
    x_bottom_children = [CoupledVariable(node, parent, -1.0) for (node, parent) in zip(model.nodes[i_right_bottom[2:end]], x_bottom_parents)]
    
    y_bottom_parents = [SpatialVariable(node, 0., dymin_bottom, dymax_bottom, :Y) for node in model.nodes[i_left_bottom]]
    y_bottom_children = [CoupledVariable(node, parent) for (node, parent) in zip(model.nodes[i_right_bottom], y_bottom_parents)]

    x_top_parents = [SpatialVariable(node, 0., dxmin, dxmax, :X) for node in model.nodes[i_left_top]]
    x_top_children = [CoupledVariable(node, parent, -1.0) for (node, parent) in zip(model.nodes[i_right_top], x_top_parents)]

    y_top_parents = [SpatialVariable(node, 0., dymin_top, dymax_top, :Y) for node in model.nodes[i_left_top]]
    y_top_children = [CoupledVariable(node, parent) for (node, parent) in zip(model.nodes[i_right_top], y_top_parents)]

    x_sup_parent = SpatialVariable(model.nodes[i_support_mid_left], 0., dxmin, dxmax, :X)
    x_sup_child = CoupledVariable(model.nodes[i_support_mid_right], x_sup_parent, -1.0)

    y_sup_parent = SpatialVariable(model.nodes[i_support_mid_left], 0., dymin_support, dymax_support, :Y)
    y_sup_child = CoupledVariable(model.nodes[i_support_mid_right], y_sup_parent)

    # bottom chord
    bottom_parent = AreaVariable(first(model.elements[:bottomchord]), .9Amax_bottom, Amin_bottom, Amax_bottom)
    bottom_children = [CoupledVariable(element, bottom_parent) for element in model.elements[:bottomchord][2:end]]

    # top chord
    top_parent = AreaVariable(first(model.elements[:topchord]), .9Amax_top, Amin_top, Amax_top)
    top_children = [CoupledVariable(element, top_parent) for element in model.elements[:topchord][2:end]]

    #bottom web
    bottomweb_parent = AreaVariable(first(model.elements[:bottomweb]), .9Amax_bottomweb, Amin_bottomweb, Amax_bottomweb)
    bottomweb_children = [CoupledVariable(element, bottomweb_parent) for element in model.elements[:bottomweb][2:end]]

    #topweb
    topweb_parent = AreaVariable(first(model.elements[:topweb]), .9Amax_topweb, Amin_topweb, Amax_topweb)
    topweb_children = [CoupledVariable(element, topweb_parent) for element in model.elements[:topweb][2:end]]

    #strut
    strut_parent = AreaVariable(first(model.elements[:strut]), .9Amax_strut, Amin_strut, Amax_strut)
    strut_children = [CoupledVariable(element, strut_parent) for element in model.elements[:strut][2:end]]

    #support
    support_parent = AreaVariable(first(model.elements[:support]), .9Amax_support, Amin_support, Amax_support)
    support_children = [CoupledVariable(element, support_parent) for element in model.elements[:support][2:end]]

    # COLLECT ALL VARIABLS
    vars = TrussVariable[
        # x_bottom_parents;
        # x_bottom_children;
        # y_bottom_parents;
        # y_bottom_children;
        # x_top_parents;
        # x_top_children;
        # y_top_parents;
        # y_top_children;
        # x_sup_parent;
        # x_sup_child;
        # y_sup_parent;
        # y_sup_child;
        bottom_parent;
        bottom_children;
        top_parent;
        top_children;
        bottomweb_parent;
        bottomweb_children;
        topweb_parent;
        topweb_children;
        strut_parent;
        strut_children;
        support_parent;
        support_children
    ]
end

# params
begin
    params = TrussOptParams(model, vars)
    x0 = copy(params.values)

    # define element-wise values
    fmin = zeros(model.nElements)
    fmax = zeros(model.nElements)
    ecc = zeros(model.nElements)
    colors = Vector{Any}(undef, model.nElements)

    for (mat, group) in zip(element_material_groups, element_id_groups)
        for i in findall(model.elements, group)

            if mat == :steel
                fmin[i] = fy_steel
                fmax[i] = fy_steel
                ecc[i] = ECC_steel
                colors[i] = :gray
            else
                fmin[i] = fc_glulam
                fmax[i] = ft_glulam
                ecc[i] = ECC_glulam
                colors[i] = kjl_orange
            end

        end
    end
end

#initial visualization
begin
    g0 = Geo(updatemodel(params, x0))
    p0 = Point2.(g0.nodes_xy)
    e0 = p0[vcat(g0.indices[i_structure]...)]

    fig = Figure()

    ax = Axis(fig[1,1], aspect = DataAspect())
    ylims!(ax, 0., 2dy + dymax_top + support_vertical_offset)

    linesegments!(
        e0,
        color = colors[i_structure]
    )

    fig
end


#=
objectives and constraints
=#

begin
    function obj(x, p)
        g = GeometricProperties(x, p)

        #individual volumes
        sum(g.A[i_structure] .* g.L[i_structure] .* ecc[i_structure])
    end

    function cstr(x, p, dmax, fc, ft)
        res = solve_truss(x, p)
        uvert = abs.(res.U[2:3:end])[i_deck_nodes]

        stresses = axial_force(res, p) ./ res.A
        
        stress_cstr = [s < 0. ? -s - flower : s - fupper for (s, flower, fupper) in zip(stresses, fc, ft)]

        return [
            (uvert .- dmax);
            maximum(stress_cstr[i_bottom_elements]);
            maximum(stress_cstr[i_top_elements]);
            maximum(stress_cstr[i_bottomweb_elements]);
            maximum(stress_cstr[i_topweb_elements]);
            maximum(stress_cstr[i_strut_elements]);
            maximum(stress_cstr[i_support_elements])
        ]
    end
end

#test
begin
    OBJ = x -> obj(x, params)
    CSTR = x -> cstr(x, params, dispmax, fmin, fmax)
    c0 = CSTR(x0)
end

all(c0 .< 0)

# optimize
begin
    alg = NLoptAlg(:LD_MMA)
    opts = NLoptOptions(
        maxtime = 120
    )

    res = constrained_optimization(params, OBJ, CSTR, alg, opts)
end

# post process
GLMakie.activate!()
begin
    @show res.obj_opt
    @show res.time
    @show n_iter = length(res)
    @show res.opt_stop_type

    x_opt = res.x_opt

    tinc = range(0, res.time, n_iter)

    gset = [Geo(updatemodel(params, x)) for x in res.x_history]
    dfac = Observable(0.)

    

    i = Observable(1)
    g = @lift(gset[$i])
    p = @lift(Point2.($g.nodes_xy .+ $dfac .* $g.disp_xy))
    e = @lift($p[ls_indices])
    lw = @lift($g.areas[i_structure] ./ maximum($g.areas[i_structure]) .* 3)

    loss = @lift(Point2.(tinc[1:$i], res.obj_history[1:$i]))
end

GLMakie.activate!()
# visualize
begin
    c = colors[i_structure]
    fig = Figure(
    )

    ax = Axis(fig[1,1], aspect = DataAspect())
    ylims!(ax, 0., 2dy + dymax_top + support_vertical_offset)

    linesegments!(els, color = (:black, .25), linewidth = .5)

    linesegments!(
        e,
        # linewidth = lw,
        color = c
    )

    axloss = Axis(
        fig[2,1],
        xlabel = "TIME [s]",
        ylabel = "EMBODIED CARBON [kg CO₂e]",
        aspect = nothing
    )

    style1!(axloss)
    
    lines!(
        loss,
        color = :black,
        label = "OBJ"
        )

    on(i) do _
        autolimits!(axloss)
    end

    rowsize!(fig.layout, 2, Aspect(1, .3))

    fig
end

i[] = 1
for k = 1:length(res)
    i[] = k
    sleep(1e-2)
end

Vwood = 0.
Vsteel = 0.

for element in res.model_opt.elements

    if element.id == :deck
        continue
    end
    if element.section.E < 100e6
        Vwood += element.section.A * element.length
    else
        Vsteel += element.section.A * element.length
    end
end

i[] = length(res)

run_code = reduce(*, [material_dict[mat] for mat in element_material_groups])
savename = date * "_" * run_code

# jldsave(datadir * savename * ".jld2"; results = res)
# GHsave(res.model_opt, datadir * savename * ".json")