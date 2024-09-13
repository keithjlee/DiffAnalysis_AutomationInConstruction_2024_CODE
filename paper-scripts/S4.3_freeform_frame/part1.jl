include("problem_formulation.jl")

#=
MAKE A MODEL THAT IS ONLY THE MAIN STRUCTURE
=#

model_init = Asap.Model(nodes[:main], e_main, loads)
Asap.solve!(model_init)

begin
    geo_init = Geo(model_init)

    p_init = Point3.(geo_init.nodes)
    e_init = p_init[geo_init.indices_flat]

    fig = Figure()
    ax = Axis3(
        fig[1,1],
        aspect = :data
    )
    hidedecorations!(ax); hidespines!(ax)

    linesegments!(e_init)
    scatter!(p_init, color = :white)

    fig
end

#=
single section sizing.

Given a tube [d, alpha], where:

d = outer diameter
alpha * d = inner diameter

Find the combination of d, alpha that minimizes structural volume while satisfying all stress constraints
=#

# setup (use of dummy params for dof indexing and K assembly)
begin
    dummy_vars = [AreaVariable(model_init.elements[1], 1e-2, 1e-4, 2e-1)]
    dummy_params = FrameOptParams(model_init, dummy_vars)

    cst_parameters = (
        nelements = model_init.nElements,
        ndofs = model_init.nDOFs,
        freeids = model_init.freeDOFs,
        L = getproperty.(model_init.elements, :length),
        E = 200e6,
        G = 80e6,
        fy = 350e3,
        R = getproperty.(model_init.elements, :R),
        params = dummy_params,
        dmax = dmax
    )
end

#objectives and constraints
begin
    function obj(x, p)

        area = a_tube(x[1], x[2])

        return reduce(+, area .* p.L)

    end


    function cstr(x, p; alg = UMFPACKFactorization())

        #section properties
        a = a_tube(x[1], x[2])
        i = I_tube(x[1], x[2])
        j = 2 * i
        s = j / x[1]

        A = fill(a, p.nelements)
        Ixy = fill(i, p.nelements)
        J = fill(j, p.nelements)

        ke = DiffAnalysis_AIC2024.k_frame.(p.E, p.G, A, p.L, Ixy, Ixy, J)

        Ke = DiffAnalysis_AIC2024.get_global_ks(p.R, ke)

        K = DiffAnalysis_AIC2024.assemble_global_K(Ke, p.params)

        u = DiffAnalysis_AIC2024.solve_u(K, p.params, alg)

        U = DiffAnalysis_AIC2024.replace_values_buffer(zeros(p.ndofs), p.freeids, u)


        #vertical displacements
        uvert = abs.(U[3:6:end])

        # stresses
        local_forces = DiffAnalysis_AIC2024.Flocal(U, Ke, p.R, p.params)
        axial_stresses = [abs(f[1]) / a for f in local_forces]
        flexural_stresses = [abs(f[6]) / s for f in local_forces]

        #peak cstr
        return [
            maximum(uvert) - p.dmax,
            maximum(axial_stresses .+ flexural_stresses) - p.fy,
            # maximum(flexural_stresses) - p.fy
        ]

    end
end

#initialize
begin
    x0 = [.75, .5]
    lb = [0.1, .05]
    ub = [1.0, 0.98]
end

#test
begin
    OBJ = x -> obj(x, cst_parameters)
    @time o0, do0 = withgradient(OBJ, x0);

    CSTR = x -> cstr(x, cst_parameters)
    @time c0, dc0 = withjacobian(CSTR, x0);
end

#check we start in feasible region
@assert all(c0 .< 0)

#=
Optimization setup
=#
begin
    alg = NLoptAlg(:LD_MMA)
    opts = NLoptOptions(
        maxeval = 1500,
        maxtime = 180
    )

    F = TraceFunction(OBJ)
    omodel = Nonconvex.Model(F)
    addvar!(omodel, lb, ub)
    add_ineq_constraint!(omodel, CSTR)
end

#optimize

begin
    t0 = time()
    res = Nonconvex.optimize(
        omodel,
        alg,
        x0,
        options = opts
    )
    restime = time() - t0
end

# results
obj_history = getproperty.(F.trace, :output)
lines(obj_history)
x_opt = res.minimizer
c_opt = CSTR(x_opt)
optimal_volume = OBJ(x_opt)

tinc = range(0, restime, length(obj_history))

using CairoMakie; CairoMakie.activate!()
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
save(figdir * "sol1_trace.pdf", fig)

#save
A_opt = a_tube(x_opt...)
I_opt = I_tube(x_opt...)
J_opt = J_tube(x_opt...)

section_opt = Section(A_opt, E, G, I_opt, I_opt, J_opt)


for element in model_init.elements
    element.section = section_opt
end

# jldsave(datadir * "08_08_2024_single_tube_section.jld2"; section = section_opt)
# GHsave(model_init, datadir * "08_08_2024_single_tube.json")
# jldsave(datadir * "08_08_2024_single_tube.jld2"; obj_history = obj_history, d_opt = x_opt[1], alpha_opt = x_opt[2])