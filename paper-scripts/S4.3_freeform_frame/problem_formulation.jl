using DiffAnalysis_AIC2024
set_theme!(aic)
using Nonconvex, NonconvexNLopt
datadir = joinpath(@__DIR__, "data/")
isdir(datadir) || mkdir(datadir)
figdir = joinpath(@__DIR__, "figures/")

data = JSON.parsefile(joinpath(@__DIR__, "08_08_2024_funicularize.json"))
#values
dmax = 50 / 300
a_tube(d, alpha) = .25pi * d^2 * (1 - alpha^2)
I_tube(d, alpha) = pi / 64 * d^4 * (1 - alpha^4)
S_tube(d, alpha) = 2 * I_tube(d, alpha) / d
J_tube(d, alpha) = 2 * I_tube(d, alpha)

begin
    nodes = [Node([x,y,z], [true, false, true, false, true, false]) for (x,y,z) in zip(data["X"], data["Y"], data["Z"])]

    # HSS 12.5 x .25 (in)
    A = 8.98 * .0254^2
    Ix = Iy = 169 * .0254^4
    J = 338 * .0254^4
    E = 200e6
    G = 80e6
    fy = 350e3

    main_sec = Section(A, E, G, Ix, Iy, J)

    # main elements
    e_main = [Element(nodes[i1+1], nodes[i2+1], main_sec, :main) for (i1, i2) in zip(data["MAIN1"], data["MAIN2"])]

    # spine elements
    Aspine = a_tube(.2, .95)
    Ispine = I_tube(.2, .95)
    Jspine = J_tube(.2, .95)

    spine_sec = Section(Aspine, E, G, Ispine, Ispine, Jspine)

    e_spine = [Element(nodes[i1+1], nodes[i2+1], spine_sec, :spine) for (i1, i2) in zip(data["SPINE1"], data["SPINE2"])]

    # strut elements
    A_strut = 4 * .0254^2
    I_strut = pi * (4 * .0254)^4 / 64
    
    strut_sec = Section(A_strut, E, G, I_strut, I_strut, 2I_strut)

    e_strut = [Element(nodes[i1+1], nodes[i2+1], strut_sec, :strut; release = :fixedfree) for (i1, i2) in zip(data["STRUT1"], data["STRUT2"])]

    # long elements
    e_long = [Element(nodes[i1+1], nodes[i2+1], main_sec, :long) for (i1, i2) in zip(data["LONG1"], data["LONG2"])]

    for node in nodes
        if node.position[3] == 0.0
            node.dof = [false, false, false, false, true, false]
        end
    end

    #node ids
    # mark the main nodes
    for i in sort(unique([data["MAIN1"]; data["MAIN2"]]) .+ 1)
        nodes[i].id = :main
    end

    for i in data["LEFTANCHORS"] .+ 1
        nodes[i].id = :leftanchor
    end

    for i in data["RIGHTANCHORS"] .+ 1
        nodes[i].id = :rightanchor
    end



    elements = [e_main; e_spine; e_strut]

    imain = unique(sort([data["MAIN1"]; data["MAIN2"]]))

    loads = [NodeForce(node, [0., 0., -40.]) for node in nodes[imain .+ 1]]

    model = Asap.Model(nodes, elements, loads)
    Asap.solve!(model)
end

