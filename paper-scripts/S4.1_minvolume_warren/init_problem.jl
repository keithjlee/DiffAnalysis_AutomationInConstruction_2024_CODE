# load environment and activate dependencies
using DiffAnalysis_AIC2024
using Nonconvex, NonconvexNLopt
set_theme!(aic)
figdir = joinpath(@__DIR__, "figures/")
datadir = joinpath(@__DIR__, "data/")

#=

direct_volume_minimization_warren_trusses.jl

This script explores how gradient-based optimization can be used to efficiently find a minimum-volume truss design solution subject to displacement and stress constraints given an initial design topology and parameterization.

We compare the speed and computationally efficiency of our method vs. a typical topology optimization based approach.

In our problem formulation, we compare the following:

1. Adjoint accelerated AD + SLSQP (Gradient based)
2. Standard AD + SLSQP (Gradient based)
3. Finite Differencing + SLSQP (Gradient based)
4. COBYLA (Gradient free)
=#

#=
PROBLEM FORMULATION
Units are in [kN], [m]
=#

# initial section
begin
    A = 0.01 #m²
    E = 200E6 #KN/m²
    fy = 350e3 #kN/m²
    section = TrussSection(A, E)
end

# Design parameters
begin
    L = 10. # total span of truss
    nx = 12 # number of bays in span
    dy = 1.0 # truss depth
    load = 75. # force applied at each free node in bottom chord (downwards)
end

# Generate truss
begin
    dx = L / nx # bay length
    dmax = L / 360 #displacement limit
    gen = Warren2D(nx, dx, dy, section; load = [0., -load, 0.]) # Warren truss generator
    model = gen.model # structural model
end

# Visualization assets
begin
    geo = Geo(model) #convert to easily plottable data structure
    dfac = Observable(0.) #this allows us to play with the displacement scale when visualizing our structures

    _nodes = @lift(Point2.(geo.nodes_xy .+ $dfac .* geo.disp_xy))
    _elements = @lift($_nodes[geo.indices_flat])
end

# Plot initial structure
begin
    fig = Figure()
    ax = Axis(
        fig[1,1],
        aspect = DataAspect(),
        xlabel = "[m]",
        ylabel = "[m]"
    )

    hidespines!(ax); tickstoggle!(ax)
    ylims!(-2dy, 2dy)

    ls_elements = linesegments!(
        _elements,
        color = :black
    )

    sc_nodes = scatter!(
        _nodes,
        color = :white
    )

    sl = Slider(
        fig[2,1],
        range = 0:0.5:100,
        startvalue = 0
    )

    on(sl.value) do val
        dfac[] = val
    end

    text!(_nodes; text = ["N$i" for i = 1:model.nNodes])
    text!(Point2.(midpoint.(model.elements)), color = kjl_blue; text = ["E$i" for i = 1:model.nElements])
    fig
end

#=
OPTIMIZATION: VARIABLES AND PARAMETERS
Defining our design variables and optimization parameters.
=#

node_pairs = find_free_node_pairs(gen)
element_pairs = find_element_pairs(gen)

# variable bounds
begin
    # Area variable parameters
    Amin = 1e-4 #m²
    Amax = 0.2 #m²
    A0 = 0.5Amax #starting guess

    # Spatial variable parameters. Note by default spatial variables are ADDITIVE.
    xmin = -0.49 * dx #m (use slightly less than half of the bay length to prevent nodes from crossing)
    xmax = 0.49 * dx #m
    x0 = 0. #starting guess

    # Top node variables
    ymin_top = -0.9 * dy
    ymax_top = dy 
    y0_top = 0. # 0.9dy

    # Bottom node variables
    ymin_bottom = -dy
    ymax_bottom = 0.75dy
    y0_bottom = .5*ymin_bottom#-0.9dy
end

# make variables (coupled)
begin
    # Design variables
    vars = [AreaVariable(element, A, Amin, Amax) for element in model.elements] # this is a weird initialization trick for segfault stability ?
    vars = TrussVariable[]

    #element areas
    push!(vars, AreaVariable(model.elements[element_pairs.i_mid], A0, Amin, Amax))
    for (il, ir) in zip(element_pairs.i_left, element_pairs.i_right)
        push!(vars, AreaVariable(model.elements[il], A0, Amin, Amax))
        push!(vars, CoupledVariable(model.elements[ir], last(vars)))
    end
    
    #node positions
    # top nodes
    for i in node_pairs.i_top_mid
        push!(vars, SpatialVariable(model.nodes[i], y0_top, ymin_top, ymax_top, :Y))
    end

    for (il, ir) in zip(node_pairs.i_top_left, node_pairs.i_top_right)

        # Y direction
        push!(vars, SpatialVariable(model.nodes[il], y0_top, ymin_top, ymax_top, :Y))
        push!(vars, CoupledVariable(model.nodes[ir], last(vars)))

        # X direction
        push!(vars, SpatialVariable(model.nodes[il], x0, xmin, xmax, :X))
        push!(vars, CoupledVariable(model.nodes[ir], last(vars), -1.0))
    end

end

# make optimization parameters and extract initial design vector
begin
    params = TrussOptParams(model, vars)
    x_init = params.values
end