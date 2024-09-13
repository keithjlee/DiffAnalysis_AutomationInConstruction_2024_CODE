#run part 2 first

model2 = res.model_opt

#variables
begin
    # main elements
    main_parent = AreaVariable(model2.elements[:frame][1], A, 1e-5, 10A)
    main_children = [CoupledVariable(element, main_parent) for element in model2.elements[:frame][2:end]]

    # spine
    spine_parent = AreaVariable(model2.elements[:spine][1], Aspine, 1e-5, 10A)
    spine_children = [CoupledVariable(element, spine_parent) for element in model2.elements[:spine][2:end]]

    # strut
    strut_parent = AreaVariable(model2.elements[:strut1][1], A_strut, 1e-5, 10A)
    strut_children = [CoupledVariable(element, strut_parent) for element in [model2.elements[:strut1][2:end]; model2.elements[:strut2]]]

    # tie
    tie_parent = AreaVariable(model2.elements[:tie][1], A_strut, 1e-5, 10A)
    tie_children = [CoupledVariable(element, tie_parent) for element in model2.elements[:tie][2:end]]

    vars2 = FrameVariable[
        main_parent;
        main_children;
        spine_parent;
        spine_children;
        strut_parent;
        strut_children;
        tie_parent;
        tie_children
    ]
end

# make parameters
begin
    params2 = FrameOptParams(model2, vars2)
    x0 = params2.values #MAIN, SPINE, STRUT

    cst_params2 = (
        L = getproperty.(model2.elements, :length),
        fy = 350e3,
        R = getproperty.(model2.elements, :R),
        dmax = 50/300,
        iframe = findall(model2.elements, :frame),
        istrut = sort([findall(model2.elements, :strut1); findall(model2.elements, :strut2)]),
        ispine = findall(model2.elements, :spine),
        itie = findall(model2.elements, :tie)
        )
end

strut_lengths = getproperty.([model2.elements[:strut1]; model2.elements[:strut2]], :length)

#initial values and bounds
begin
    d0_main = d_init #m
    a0_main = 0.5

    d0_spine = d_init #m
    a0_spine = 0.5

    d0_strut = d_init #m
    a0_strut = 0.5

    d0_tie = d_init
    a0_tie = 0.5

    amin = 0.01
    amax = 0.98

    dmin = 0.2 #m
    dmax = 1.0 #m

    x0 = [d0_main, a0_main, d0_spine, a0_spine, d0_strut, a0_strut, d0_tie, a0_tie]
    # lb = [dmin1, amin, dmin1, amin, dmin2, amin, dmin2, amin]
    # ub = [dmax1, amax, dmax1, amax, dmax2, amax, dmax2, amax]
    lb = repeat([dmin, amin], 4)
    ub = repeat([dmax, amax], 4)
end

# objective, constraint, and closures
begin
    function cstr(x, p, p2; alg = UMFPACKFactorization())

        # make section propertues
        dmain, alphamain, dspine, alphaspine, dstrut, alphastrut, dtie, alphatie = x

        a_main = a_tube(dmain, alphamain)
        i_main = I_tube(dmain, alphamain)
        j_main = 2i_main
        s_main = j_main / dmain

        a_spine = a_tube(dspine, alphaspine)
        i_spine = I_tube(dspine, alphaspine)
        j_spine = 2i_spine
        s_spine = j_spine / dspine

        a_strut = a_tube(dstrut, alphastrut)
        i_strut = I_tube(dstrut, alphastrut)
        j_strut = 2i_strut
        s_strut = j_strut / dstrut

        a_tie = a_tube(dtie, alphatie)
        i_tie = I_tube(dtie, alphatie)
        j_tie = 2i_tie
        s_tie = j_tie / dtie

        newA = [a_main, a_spine, a_strut, a_tie]
        newI = [i_main, i_spine, i_strut, i_tie]
        newJ = [j_main, j_spine, j_strut, j_tie]
        newS = [s_main, s_spine, s_strut, s_tie]

        A = DiffAnalysis_AIC2024.replace_values(p.A, p.indexer.iA, newA[p.indexer.iAg])
        Ixy = DiffAnalysis_AIC2024.replace_values(p.Ix, p.indexer.iA, newI[p.indexer.iAg])
        J = DiffAnalysis_AIC2024.replace_values(p.J, p.indexer.iA, newJ[p.indexer.iAg])
        S = newS[p.indexer.iAg]

        ke = DiffAnalysis_AIC2024.k_frame.(p.E, p.G, A, p2.L, Ixy, Ixy, J)

        Ke = DiffAnalysis_AIC2024.get_global_ks(p2.R, ke)

        K = DiffAnalysis_AIC2024.assemble_global_K(Ke, p)

        u = DiffAnalysis_AIC2024.solve_u(K, p, alg)

        U = DiffAnalysis_AIC2024.replace_values(zeros(p.n), p.freeids, u)

        #vertical displacements
        uvert = abs.(U[3:6:end])

        # stresses
        local_forces = DiffAnalysis_AIC2024.Flocal(U, Ke, p2.R, p)
        stresses = abs.(getindex.(local_forces, 1)) ./ A .+ abs.(getindex.(local_forces, 6)) ./ S

        return [
            maximum(uvert) - p2.dmax,
            maximum(stresses[p2.iframe]) - p2.fy,
            maximum(stresses[p2.istrut]) - p2.fy,
            maximum(stresses[p2.ispine]) - p2.fy,
            maximum(stresses[p2.itie]) - p2.fy,
            dstrut - dspine,
            dstrut - dmain,
            dtie - dspine
        ]
    end

    function obj(x, p, p2)

        # make section propertues
        dmain, alphamain, dspine, alphaspine, dstrut, alphastrut, dtie, alphatie = x
        
        a_main = a_tube(dmain, alphamain)
        a_spine = a_tube(dspine, alphaspine)
        a_strut = a_tube(dstrut, alphastrut)
        a_tie = a_tube(dtie, alphatie)
        newA = [a_main, a_spine, a_strut, a_tie]

        A = DiffAnalysis_AIC2024.replace_values(p.A, p.indexer.iA, newA[p.indexer.iAg])

        dot(A, p2.L)
    end

    OBJ = x -> obj(x, params2, cst_params2)
    o0, do0 = withgradient(OBJ, x0)
    CSTR = x -> cstr(x, params2, cst_params2)
    @time c0, dc0 = withjacobian(CSTR, x0)
end

c0
@assert all(c0 .<= 0)

# opt params
begin
    alg = NLoptAlg(:LD_MMA)
    opts = NLoptOptions(
        maxeval = 1000,
        maxtime = 300
    )
end

# solve
begin
    F = TraceFunction(OBJ)
    omodel = Nonconvex.Model(F)
    addvar!(omodel, lb, ub)
    add_ineq_constraint!(omodel, CSTR)

    t0 = time()
    res = Nonconvex.optimize(
        omodel,
        alg,
        x0,
        options = opts
    )
    restime = time() - t0
end

@show res.minimizer
CSTR(res.minimizer)

#post process
obj_history = getproperty.(F.trace, :output)
tinc = range(0, restime, length(obj_history))

lines(tinc, obj_history)

opt_model = updatemodel(params2, res.minimizer)

begin
    fig = Figure(size = halfwidth())
    ax = Axis(
        fig[1,1],
        xlabel = "TIME [s]",
        ylabel = "VOLUME [m³]",
        aspect = 2.5
    )

    style1!(ax)

    lines!(tinc, obj_history, color = kjl_blue)

    fig
end

dframe, aframe, dspine, aspine, dstrut, astrut, dtie, atie = res.minimizer

# jldsave(datadir * "08_12_2024_final2.jld2"; model = opt_model, dframe = dframe, aframe = aframe, dspine = dspine, aspine = aspine, dstrut = dstrut, astrut = astrut, dtie = dtie, atie = atie, obj_history = obj_history, time_history = tinc)