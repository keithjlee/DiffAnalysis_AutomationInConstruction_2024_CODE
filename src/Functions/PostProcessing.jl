"""
    Flocal(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)
[2×1] vector of end element end forces in LCS
"""
function Flocal(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::AbstractOptParams)

    #Displacements w/r/t each element
    uEs = [u[id] for id in p.dofids]

    #End forces
    Rs .* Eks .* uEs
end

function ChainRulesCore.rrule(::typeof(Flocal), u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)

    F = Flocal(u, Eks, Rs, p)

    function Flocal_pullback(F̄)
        du = zero(u)
        dK = zero.(Eks)
        dR = zero.(Rs)
        ids = p.dofids

        for i in eachindex(ids)
            # F̄ ⋅ dF/dR
            dR[i] = F̄[i] * (Eks[i] * u[ids[i]])'

            #F̄ ⋅ dF/du
            du[ids[i]] .+= (Rs[i] * Eks[i])' * F̄[i]

            #F̄ ⋅ dF/dK
            dK[i] = (Rs[i]' * F̄[i]) * u[ids[i]]'
        end

        return (NoTangent(), du, dK, dR, NoTangent())

    end

    return F, Flocal_pullback
end

function Flocal_noadjoint(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::AbstractOptParams)

    #Displacements w/r/t each element
    uEs = [u[id] for id in p.dofids]

    #End forces
    Rs .* Eks .* uEs
end

function element_forces(res::FrameResults, p::FrameOptParams)
    return Flocal(res.U, res.K, res.R, p)
end

"""
    Flocal(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)
[2×1] vector of end element end forces in LCS
"""
function axial_force(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)
    fvecs = Flocal(u, Eks, Rs, p)
    getindex.(fvecs, 2)
end

function axial_force_noadjoint(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)
    fvecs = Flocal_noadjoint(u, Eks, Rs, p)
    getindex.(fvecs, 2)
end

function axial_force(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::FrameOptParams)
    fvecs = Flocal(u, Eks, Rs, p)
    getindex.(fvecs, 7)
end

"""
    σaxial(F::Vector{Float64}, A::Vector{Float64})

Get the absolute value of the axial stress experienced under force F and area A
"""
function σaxial(F::Vector{Float64}, A::Vector{Float64})
    return abs.(F) ./ A
end

"""
dg/dF = dg/dσ ⋅ dσ/dF = σ̄ ⋅ F / |F| / A
dg/dA = dg/dσ ⋅ dσ/dA = σ̄ ⋅ -F / A²
"""
function ChainRulesCore.rrule(::typeof(σaxial), F::Vector{Float64}, A::Vector{Float64})
    stress = σaxial(F, A)

    function σaxial_pullback(σ̄)
        s = sign.(F)
        return (NoTangent(), σ̄ .* s ./ A, -σ̄ .* s .* F ./ A.^2)
    end

    return stress, σaxial_pullback
end

"""
    axialforce(t::TrussResults, p::TrussOptParams)
Axial forces in truss structure
"""
function axial_force(t::TrussResults, p::TrussOptParams)
    axial_force(t.U, t.K, t.R, p)
end

"""
    axial_stress(t::TrussResults, p::TrussOptParams)
Axial stresses in truss structure
"""
function axial_stress(t::TrussResults, p::TrussOptParams)
    F = axial_force(t.U, t.K, t.R, p)
    σaxial(F, t.A)
end

"""
    axialforce(t::NetworkResults, p::NetworkOptParams)
Axial forces in FDM network
"""
function axial_force(t::NetworkResults, p::NetworkOptParams)
    norm.(eachrow(p.C * [t.X t.Y t.Z])) .* t.Q
end

export f_axial_new
function f_axial_new(U::Vector{Float64}, Ks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)

    # initialize
    forces = Vector{Float64}(undef, length(Ks))

    for i in eachindex(forces)
        edof = p.dofids[i]
        edisp = U[edof]

        forces[i] = dot(Rs[i][2,:], Ks[i] * edisp)
    end

    return forces

end

function ChainRulesCore.rrule(::typeof(f_axial_new), U::Vector{Float64}, Ks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)

    F = f_axial_new(U, Ks, Rs, p)

    function f_axial_new_pullback(F̄)

        dU = zero(U)
        dKs = zero.(Ks)
        dRs = zero.(Rs)

        for i in eachindex(Ks)

            r = Rs[i][2, :]
            k = Ks[i]
            u = U[p.dofids[i]]

            dU[p.dofids[i]] .+= F̄[i] .* k * r
            dKs[i] = F̄[i] .* r * u'
            dRs[i][2, :] = F̄[i] .* u'* k

            return (NoTangent(), dU, dKs, dRs, NoTangent())
        end

    end

    return F, f_axial_new_pullback
end
"""
    updatemodel(p::TrussOptParams, u::Vector{Float64})

Generate a new structural model from the results of an optimization
"""
function updatemodel(p::TrussOptParams, u::Vector{Float64})
    
    #final values
    X = add_values(p.X, p.indexer.iX, u[p.indexer.iXg] .* p.indexer.fX)
    Y = add_values(p.Y, p.indexer.iY, u[p.indexer.iYg] .* p.indexer.fY)
    Z = add_values(p.Z, p.indexer.iZ, u[p.indexer.iZg] .* p.indexer.fZ)
    A = replace_values(p.A, p.indexer.iA, u[p.indexer.iAg] .* p.indexer.fA)

    #new model
    nodes = Vector{TrussNode}()
    elements = Vector{TrussElement}()
    loads = Vector{NodeForce}()

    #new nodes
    for (node, x, y, z) in zip(p.model.nodes, X, Y, Z)
        newnode = TrussNode([x, y, z], node.dof)
        newnode.id = node.id
        push!(nodes, newnode)
    end

    #new elements
    for (id, e, a, el) in zip(p.nodeids, p.E, A, p.model.elements)
        newelement = TrussElement(nodes[id]..., TrussSection(a, e))
        newelement.id = el.id
        push!(elements, newelement)
    end

    #new loads
    for load in p.model.loads
        newload = NodeForce(nodes[load.node.nodeID], load.value)
        newload.id = load.id
        push!(loads, newload)
    end

    model = TrussModel(nodes, elements, loads)
    Asap.solve!(model)

    return model

end

"""
    updatenetwork(p::TrussOptParams, u::Vector{Float64})

Generate a new FDM network from the results of an optimization
"""
function updatenetwork(p::NetworkOptParams, u::Vector{Float64})
    
    #final values
    X = add_values(p.X, p.indexer.iX, u[p.indexer.iXg] .* p.indexer.fX)
    Y = add_values(p.Y, p.indexer.iY, u[p.indexer.iYg] .* p.indexer.fY)
    Z = add_values(p.Z, p.indexer.iZ, u[p.indexer.iZg] .* p.indexer.fZ)
    Q = replace_values(p.q, p.indexer.iQ, u[p.indexer.iQg] .* p.indexer.fQ)

    #new model
    nodes = Vector{FDMnode}()
    elements = Vector{FDMelement}()
    loads = Vector{FDMload}()

    #new nodes
    for (node, x, y, z) in zip(p.network.nodes, X, Y, Z)
        newnode = FDMnode([x, y, z], node.dof)
        newnode.id = node.id
        push!(nodes, newnode)
    end

    #new elements
    for (q, el) in zip(Q, p.network.elements)
        newelement = FDMelement(nodes, el.iStart, el.iEnd, q)
        newelement.id = el.id
        push!(elements, newelement)
    end

    #new loads
    for load in p.network.loads
        newload = FDMload(nodes[load.point.nodeID], load.force)
        push!(loads, newload)
    end

    network = Network(nodes, elements, loads)
    Asap.solve!(network)

    return network

end

release_type_to_symbol(::Type{Asap.FixedFixed}) = :fixedfixed
release_type_to_symbol(::Type{Asap.FixedFree}) = :fixedfree
release_type_to_symbol(::Type{Asap.FreeFixed}) = :freefixed
release_type_to_symbol(::Type{Asap.Joist}) = :joist

"""
    updatemodel(p::TrussOptParams, u::Vector{Float64})

Generate a new structural model from the results of an optimization
"""
function updatemodel(p::FrameOptParams, u::Vector{Float64})
    
    #final values
    X = add_values(p.X, p.indexer.iX, u[p.indexer.iXg] .* p.indexer.fX)
    Y = add_values(p.Y, p.indexer.iY, u[p.indexer.iYg] .* p.indexer.fY)
    Z = add_values(p.Z, p.indexer.iZ, u[p.indexer.iZg] .* p.indexer.fZ)

    A = replace_values(p.A, p.indexer.iA, u[p.indexer.iAg] .* p.indexer.fA)
    Ix = replace_values(p.Ix, p.indexer.iIx, u[p.indexer.iIxg] .* p.indexer.fIx)
    Iy = replace_values(p.Iy, p.indexer.iIy, u[p.indexer.iIyg] .* p.indexer.fIy)
    J = replace_values(p.J, p.indexer.iJ, u[p.indexer.iJg] .* p.indexer.fJ)

    #new model
    nodes = Vector{Node}()
    elements = Vector{Element}()
    loads = Vector{Asap.AbstractLoad}()

    #new nodes
    for (node, x, y, z) in zip(p.model.nodes, X, Y, Z)
        newnode = Node([x, y, z], node.dof, node.id)
        push!(nodes, newnode)
    end

    #new elements
    for (id, el, a, ix, iy, j, r) in zip(p.nodeids, p.model.elements, A, Ix, Iy, J, p.releases)
        section = deepcopy(el.section)

        newsection = Section(a, section.E, section.G, ix, iy, j)
        newelement = Element(nodes[id]..., newsection; release = release_type_to_symbol(r))
        newelement.id = el.id
        push!(elements, newelement)
    end

    #new loads
    for load in p.model.loads

        if typeof(load) <: Asap.NodeLoad
            newload = typeof(load)(nodes[load.node.nodeID], load.value)
            newload.id = load.id
            push!(loads, newload)
        else
            if typeof(load) == Asap.PointLoad
                newload = PointLoad(elements[load.element.elementID], load.position, load.value)
                newload.id = load.id
                push!(loads, newload)
            else
                newload = typeof(load)(elements[load.element.elementID], load.value)
                newload.id = load.id
                push!(loads, newload)
            end
        end
    end

    model = Model(nodes, elements, loads)
    Asap.solve!(model)

    return model

end

export FL2
function FL2(F, L)

    loss = 0.0
    for (f, l) in zip(F, L)
        if f ≥ 0
            loss += f * l
        else
            loss -= f * l^2
        end
    end

    return loss
end

function ChainRulesCore.rrule(::typeof(FL2), F, L)

    loss = FL2(F, L)

    function FL2_pullback(Ō)

        F̄ = zero(F)
        L̄ = zero(L)

        i_tens = findall(F .≥ 0)
        i_comp = findall(F .< 0)

        F̄[i_tens] .+= L̄[i_tens] .* Ō
        F̄[i_comp] .-=  (L̄[i_comp]).^2 .* Ō

        L̄[i_tens] .+= F̄[i_tens] .* Ō
        L̄[i_comp] .-= 2F̄[i_comp] .* L̄[i_comp] .* Ō

        return NoTangent(), F̄, L̄
    end

    return loss, FL2_pullback
end

export FL2_noadjoint
function FL2_noadjoint(F, L)

    loss = 0.0
    for (f, l) in zip(F, L)
        if f ≥ 0
            loss += f * l
        else
            loss -= f * l^2
        end
    end

    return loss
end