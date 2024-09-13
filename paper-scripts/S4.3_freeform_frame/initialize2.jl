using DiffAnalysis_AIC2024
set_theme!(aic)
using Nonconvex, NonconvexNLopt
datadir = joinpath(@__DIR__, "data/")
isdir(datadir) || mkdir(datadir)
figdir = joinpath(@__DIR__, "figures/")

#=
READ AND GENERATE STRUCTURE
=#

begin

    E = 200e6
    G = 80e6

    data = JSON.parsefile(joinpath(@__DIR__, "08_12_2024_funicularize2.json"))

    I1 = data["I1"] .+ 1
    I2 = data["I2"] .+ 1

    #values
    dmax = 50 / 300
    a_tube(d, alpha) = .25pi * d^2 * (1 - alpha^2)
    I_tube(d, alpha) = pi / 64 * d^4 * (1 - alpha^4)
    S_tube(d, alpha) = 2 * I_tube(d, alpha) / d
    J_tube(d, alpha) = 2 * I_tube(d, alpha)

    nodes = [Node([x,y,z], [true, false, true, false, true, false]) for (x,y,z) in zip(data["X"], data["Y"], data["Z"])]

    # main frame elements
    d_init = .2
    alpha_init = .95
    Ainit = a_tube(d_init, alpha_init)
    Iinit = I_tube(d_init, alpha_init)
    Jinit = J_tube(d_init, alpha_init)

    sec_init = Section(Ainit, E, G, Iinit, Iinit, Jinit)
    
    e_main = [Element(nodes[I1[i]], nodes[I2[i]], sec_init, :frame) for i in data["IFRAME"] .+ 1]

    # spine elements
    e_spine = [Element(nodes[I1[i]], nodes[I2[i]], sec_init, :spine) for i in data["ISPINE"] .+ 1]

    # strut elements
    e_strut1 = [Element(nodes[I1[i]], nodes[I2[i]], sec_init, :strut1) for i in data["ISTRUT1"] .+ 1]
    e_strut2 = [Element(nodes[I1[i]], nodes[I2[i]], sec_init, :strut2) for i in data["ISTRUT2"] .+ 1]

    # tie elements
    e_tie = [Element(nodes[I1[i]], nodes[I2[i]], sec_init, :tie) for i in data["ITIE"] .+ 1]

    #support nodes
    for i in data["ISUPPORT"] .+ 1
        nodes[i].dof = [false, false, false, false, true, false]
        nodes[i].id = :support
    end

    #anchor nodes
    for i in data["IANCHORS"] .+ 1
        nodes[i].dof = [false, false, false, false, true, false]

        if nodes[i].position[1] < 0
            nodes[i].id = :leftanchor
        else
            nodes[i].id = :rightanchor
        end
    end

    #active node
    for i in I2[data["ISTRUT1"].+ 1]
        nodes[i].id = :active1
        nodes[i].dof = [true, true, true, true, true, true]
    end

    for i in I2[data["ISTRUT2"].+ 1]
        nodes[i].id = :active2
        nodes[i].dof = [true, true, true, true, true, true]
    end

    elements = [e_main; e_spine; e_strut1; e_strut2; e_tie]
    loads = [NodeForce(node, [0., 0., -40.]) for node in nodes[data["ILOAD"] .+ 1]]
    model = Asap.Model(nodes, elements, loads)
    Asap.solve!(model)
end


# visualize
begin
    geo = Geo(model)
    dfac = Observable(0.)
    pts = @lift(Point3.(geo.nodes .+ $dfac .* geo.disp))
    els = @lift($pts[geo.indices_flat])

    fig = Figure()
    ax = LScene(
        fig[1,1],
        show_axis = false
    )

    linesegments!(els)

    fig
end