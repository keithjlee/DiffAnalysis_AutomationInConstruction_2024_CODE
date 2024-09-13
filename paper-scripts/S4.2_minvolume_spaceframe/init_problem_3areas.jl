using DiffAnalysis_AIC2024
set_theme!(aic)
using Nonconvex, NonconvexNLopt

datadir = joinpath(@__DIR__, "data/"); isdir(datadir) || mkdir(datadir)
figdir = joinpath(@__DIR__, "figures/"); isdir(figdir) || mkdir(figdir)


#=
INITIALIZE STRUCTURE
=#
n = 8
d = 3.
dz = 2.25

A0 = 1e-2
E = 200e6
sec = TrussSection(A0, E)
fy = 350e3

p = 30.
load = [0., 0., -p]

ddiv = 300
dmax = n * d / ddiv

gen = SpaceFrame(n, d, n, d, dz, sec; load = load, support = :x)
model = gen.model

for node in model.nodes[:support]
    fixnode!(node, :free)
    node.id = :bottom
end

for i in [gen.ibottom[:, 1]; gen.ibottom[1,:]]
    fixnode!(model.nodes[i], :pinned)
    model.nodes[i].id = :support
end

for i in gen.itop[:, 1]
    model.nodes[i].position .-= [0., .5d, 0.]
    fixnode!(model.nodes[i], :pinned)
    model.nodes[i].id = :support
end

for i in gen.itop[1, :]
    model.nodes[i].position .-= [.5d, 0., 0.]
    fixnode!(model.nodes[i], :pinned)
    model.nodes[i].id = :support
end

Asap.solve!(model; reprocess = true)

newloads = [NodeForce(node, load) for node in model.nodes[:bottom]]
Asap.solve!(model, newloads)

# GHsave(model, datadir * "init.json")
geo = Geo(model)

dfac = Observable(0.)
_nodes = @lift(Point3.(geo.nodes .+ $dfac .* geo.disp))
_elements = @lift($_nodes[geo.indices_flat])


#=
OPTIMIZATION VARIABLES AND PARAMETERS
=#


grid = gen.itop
@assert iseven(n)

nmid = Int(n / 2)

iparent = grid[1:nmid, 1:nmid]
ichild1 = reverse(grid[nmid+1:end, 1:nmid], dims = 1)
fac1 = [-1., 1.]
ichild2 = reverse(grid[1:nmid, nmid+1:end], dims = 2)
fac2 = [1., -1.]
ichild3 = reverse(grid[nmid+1:end, nmid+1:end])
fac3 = [-1., 1.]

vars = [SpatialVariable(node, 0. -.5*dz, 2dz, :Z) for node in model.nodes[:top]]

vars = Vector{TrussVariable}()

zmintop = -0.9dz
zmaxtop = dz

Amin = 1e-3
Amax = 0.2

# independent variable diagonals
for i in diag(gen.itop)
    if model.nodes[i].id == :support
        continue
    end

    push!(vars, SpatialVariable(model.nodes[i], 0., zmintop, zmaxtop, :Z))
end



# top variables coupled
for i = 2:size(gen.itop, 1)
    for j = 1:i-1

        iparent = gen.itop[i,j]
        ichild = gen.itop[j,i]

        if model.nodes[iparent].id == :support || model.nodes[ichild].id == :support
            continue
        end

        parent = SpatialVariable(model.nodes[iparent], 0., zmintop, zmaxtop, :Z)
        child = CoupledVariable(model.nodes[ichild], parent)

        push!(vars, parent, child)

    end
end


itop = findall(model.elements, :top)
iweb = findall(model.elements, :web)
ibottom = findall(model.elements, :bottom)

top_parent = AreaVariable(model.elements[first(itop)], 0.5Amax, Amin, Amax)
top_children = [CoupledVariable(model.elements[i], top_parent) for i in itop[2:end]]
push!(vars, top_parent, top_children...)

web_parent = AreaVariable(model.elements[first(iweb)], 0.5Amax, Amin, Amax)
web_children = [CoupledVariable(model.elements[i], web_parent) for i in iweb[2:end]]
push!(vars, web_parent, web_children...)

bottom_parent = AreaVariable(model.elements[first(ibottom)], 0.5Amax, Amin, Amax)
bottom_children = [CoupledVariable(model.elements[i], bottom_parent) for i in ibottom[2:end]]
push!(vars, bottom_parent, bottom_children...)


i_stressed_elements = Vector{Int64}()
for i in eachindex(model.elements)
    if !(model.elements[i].nodeStart.id == model.elements[i].nodeEnd.id == :support)
        push!(i_stressed_elements, i)
    end
end

params = TrussOptParams(model, vars)
x_init = params.values