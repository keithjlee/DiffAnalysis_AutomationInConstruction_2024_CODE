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
    A = 8.98 * .0254^2
    Ix = Iy = 169 * .0254^4
    J = 338 * .0254^4
    E = 200e6
    G = 80e6
    fy = 350e3

    main_sec = Section(A, E, G, Ix, Iy, J)
    e_main = [Element(nodes[I1[i]], nodes[I2[i]], main_sec, :frame) for i in data["IFRAME"] .+ 1]

    # spine elements
    Aspine = a_tube(.2, .95)
    Ispine = I_tube(.2, .95)
    Jspine = J_tube(.2, .95)

    spine_sec = Section(Aspine, E, G, Ispine, Ispine, Jspine)
    e_spine = [Element(nodes[I1[i]], nodes[I2[i]], spine_sec, :spine) for i in data["ISPINE"] .+ 1]

    # strut elements
    d_strut = 4 * .0254
    A_strut = pi * d_strut^2 / 4
    I_strut = pi * (d_strut)^4 / 64

    strut_sec = Section(A_strut, E, G, I_strut, I_strut, 2I_strut)
    e_strut1 = [Element(nodes[I1[i]], nodes[I2[i]], strut_sec, :strut1) for i in data["ISTRUT1"] .+ 1]
    e_strut2 = [Element(nodes[I1[i]], nodes[I2[i]], strut_sec, :strut2) for i in data["ISTRUT2"] .+ 1]

    # tie elements
    d_tie = 2 * .0254
    A_tie = pi * d_tie^2 / 4
    I_tie = pi * (d_tie)^4 / 64

    tie_sec = Section(A_tie, E, G, I_tie, I_tie, 2I_tie)
    e_tie = [Element(nodes[I1[i]], nodes[I2[i]], tie_sec, :tie) for i in data["ITIE"] .+ 1]

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