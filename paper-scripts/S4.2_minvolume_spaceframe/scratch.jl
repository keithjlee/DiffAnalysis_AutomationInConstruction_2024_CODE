using DiffAnalysis_AIC2024
set_theme!(aic)

data = jldopen(joinpath(@__DIR__, "data/08_07_2024_gb_mma_3areas.jld2"))
res = data["results"]

model = res.model_opt

eids = getproperty.(model.elements, :id)
areas = getproperty.(getproperty.(model.elements, :section), :A)

unique_areas = unique(areas)