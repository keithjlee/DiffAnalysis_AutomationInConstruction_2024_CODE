include("init.jl")
date = "08_20_2024"

using CairoMakie; CairoMakie.activate!()

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
    dymax_top = 10.

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

# optimization params
alg = NLoptAlg(:LD_MMA)
opts = NLoptOptions(
    maxtime = 120
)


bottomchord_material = :steel
topchord_material = :steel
bottomweb_material = :steel
topweb_material = :steel
strut_material = :steel
support_material = :steel


#re-initialize structure
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

#apply bounds
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

# make parameters
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

#OBJECTIVES AND CONSTRAINTS
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

#CLOSURES
OBJ = x -> obj(x, params)
CSTR = x -> cstr(x, params, dispmax, fmin, fmax)

#SOLVE
res = constrained_optimization(params, OBJ, CSTR, alg, opts)

@show res.obj_opt
@show res.time
@show n_iter = length(res)
@show res.opt_stop_type

#DRAW
c = colors[i_structure]


geo = Geo(res.model_opt)
pts = Point2.(geo.nodes_xy)
els = pts[ls_indices]

fig = Figure(
    
)
ax = Axis(fig[1,1], aspect = DataAspect())

hidedecorations!(ax); hidespines!(ax)
ylims!(ax, -dy, 2dy + dymax_top + support_vertical_offset)

linesegments!(static_deck_points, color = :black, linewidth = 1)

linesegments!(
    els,
    color = c,
    linewidth = 1.25
)

scatter!(pts, color = :white, markersize = 3, strokewidth = .5)

#SAVE
run_code = reduce(*, [material_dict[mat] for mat in element_material_groups])

println("$run_code COMPLETE; SAVING")

savename = date * "_" * run_code * "_sizeonly"

jldsave(datadir * savename * ".jld2"; results = res)
GHsave(res.model_opt, datadir * savename * ".json")
save(figdir * savename * ".pdf", fig)

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

push!(codes, run_code)
push!(EC_values, res.obj_opt)
push!(vol_steel, Vsteel)
push!(vol_wood, Vwood)
