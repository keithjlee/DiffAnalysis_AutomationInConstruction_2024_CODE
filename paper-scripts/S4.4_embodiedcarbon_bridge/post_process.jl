include("init.jl")
date = "08_16_2024"
using CairoMakie
CairoMakie.activate!()

abbreviations = ["S", "W"]
code_combinations = Vector{String}()
for a in abbreviations
    for b in abbreviations
        for c in abbreviations
            for d in abbreviations
                for e in abbreviations
                    for f in abbreviations
                        push!(code_combinations, a * b * c * d * e * f)
                    end
                end
            end
        end
    end
end

file_names = [datadir * date * "_" * code * ".jld2" for code in code_combinations]
truss_results = [jldopen(file)["results"] for file in file_names]

ec_direct = Vector{Float64}()
mass_direct = Vector{Float64}()
compliance_direct = Vector{Float64}()

for res in truss_results
    Vsteel = 0.
    Vwood = 0.

    for element in res.model_opt.elements
        if element.id == :deck
            continue
        end

        if element.section.E < 199e6
            Vwood += element.section.A * element.length
        else
            Vsteel += element.section.A * element.length
        end
    end

    push!(ec_direct, Vsteel * ECC_steel + Vwood * ECC_glulam)
    push!(compliance_direct, res.model_opt.compliance)

    push!(mass_direct, Vsteel * ρ_steel + Vwood * ρ_glulam)
end

total_time = sum(getproperty.(truss_results, :time))

CairoMakie.activate!()
begin
    ms = 5
    fig = Figure(size = halfwidth(.5))

    ax1 = Axis(fig[1,1], xlabel = "COMPLIANCE [kNm]", ylabel = "EC [tonnes CO₂e]")
    style1!(ax1)
    scatter!(compliance_direct, ec_direct ./ 1e3, markersize = ms, color = :black)
    xlims!(0, 450)
    ylims!(low = 0)

    ax1.xticks = 0:100:450

    ax2 = Axis(fig[1,2], xlabel = "MASS [tonnes]")
    style1!(ax2)
    hideydecorations!(ax2)
    scatter!(mass_direct ./ 1e3, ec_direct ./ 1e3, markersize = ms, color = :black)
    xlims!(low = 0)
    ylims!(low = 0)

    fig
end

save(figdir * "mass_ec_compliance_zeroed.pdf", fig)

code_combinations[argmin(mass_direct)]
names_by_performance = code_combinations[sortperm(ec_direct)]

figdir = joinpath(@__DIR__, "processed_figures/")
isdir(figdir) || mkdir(figdir)

for i in eachindex(code_combinations)
    res = truss_results[i]
    n = code_combinations[i]
    ec = Int(round(ec_direct[i], digits = 0))

    c = [e.section.E < 100e6 ? kjl_orange : :gray for e in res.model_opt.elements if e.id != :deck]

    geo = Geo(res.model_opt)
    pts = Point2.(geo.nodes_xy)
    els = pts[ls_indices]

    fig = Figure(backgroundcolor = :transparent)
    ax = Axis(fig[1,1], aspect = DataAspect())
    ax.titlesize = 28

    hidedecorations!(ax); hidespines!(ax)
 
    # ax.titlevisible = true
    ylims!(ax, -dy, 2dy + 10 + support_vertical_offset)

    linesegments!(static_deck_points, color = :black, linewidth = 1)

    linesegments!(
        els,
        color = c,
        linewidth = 1.25
    )

    scatter!(pts, color = :white, markersize = 3, strokewidth = .5)

    fig
    # save(figdir * date * "_" * n * "_notitle.pdf", fig)
end

i_global_optimal = argmin(ec_direct)

begin
    fig = Figure(size = halfwidth(; factor = 1.))

    ax1 = Axis(fig[1,1], xlabel = "COMPLIANCE [kNm]", ylabel = "EMBODIED CARBON [kg CO₂e]", aspect = 2)
    style1!(ax1)
    scatter!(compliance_direct, ec_direct, markersize = 2, color = :black)
    labelscale!(ax1, 0.5)

    ax2 = Axis(fig[2,1], xlabel = "MASS [kg]", ylabel = "EMBODIED CARBON [kg CO₂e]", aspect = 2)
    style1!(ax2)
    scatter!(mass_direct, ec_direct, markersize = 2, color = :black)

    labelscale!(ax2, 0.5)

    fig
end

# save(figdir * "mass_ec_compliance.pdf", fig)


# split indices into wood and steel
i_wood = Vector{Vector{Int64}}()
i_steel = Vector{Int64}()

testres = first(truss_results)
for i in eachindex(testres.model_opt.elements)
    e = testres.model_opt.elements[i]

    if e.section.E < 199e6
        push!(i_wood, i)
    else
        push!(i_steel, i)
    end
end