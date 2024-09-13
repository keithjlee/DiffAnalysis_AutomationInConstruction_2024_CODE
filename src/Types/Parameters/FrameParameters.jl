function get_element_release(element::Element{R}) where {R}
    return R
end

struct FrameOptParams{V} <: AbstractOptParams
    model::Model #the reference  model for optimization
    values::Vector{Float64} #design variables
    indexer::FrameOptIndexer #pointers to design variables and full variables
    variables::Vector{V}
    releases
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    Ψ::Vector{Float64} #all roll angles |n_element|
    E::Vector{Float64} #all element young's modulii |n_element|
    G::Vector{Float64} #all element shear modulii |n_element|
    A::Vector{Float64} #all element areas |n_element|
    Ix::Vector{Float64} #all element X-X moment of inertia (nominal strong) |n_element|
    Iy::Vector{Float64} #all element Y-Y moment of inertia (nominal weak) |n_element|
    J::Vector{Float64} #all element torsion constant |n_element|
    P::Vector{Float64} # External load vector
    Pf::Vector{Float64} # Fixed end force vector
    C::SparseMatrixCSC{Int64, Int64} #connectivity matrix
    lb::Vector{Float64} #lower bounds of variables
    ub::Vector{Float64} #upper bounds of variables
    cp::Vector{Int64} #S.colptr
    rv::Vector{Int64} #S.rowval
    nnz::Int64 #length(S.nzval)
    inzs::Vector{Vector{Int64}} # Indices of elemental K in global S.nzval
    freeids::Vector{Int64} # [DofFree1, DofFree2,...]
    nodeids::Vector{Vector{Int64}} # [[iNodeStart, iNodeEnd] for element in elements]
    dofids::Vector{Vector{Int64}} # [[dofStartNode..., dofEndNode...] for element in elements]
    n::Int64 #total number of DOFs

    function FrameOptParams(model::Model, variables::Vector{T}) where T <: FrameVariable

        #assert model is processed
        model.processed || (Asap.process!(model))

        #global parameters
        
        #spatial positions
        xyz = Matrix{Float64}(undef, model.nNodes, 3)
        for (i, node) in enumerate(model.nodes)
            xyz[i, :] .= node.position
        end

        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]

        #element-wise collectors
        Ψ = Vector{Float64}(undef, model.nElements)
        E = Vector{Float64}(undef, model.nElements)
        G = Vector{Float64}(undef, model.nElements)
        A = Vector{Float64}(undef, model.nElements)
        Ix = Vector{Float64}(undef, model.nElements)
        Iy = Vector{Float64}(undef, model.nElements)
        J = Vector{Float64}(undef, model.nElements)
        nodeids = Vector{Vector{Int64}}(undef, model.nElements)
        dofids = Vector{Vector{Int64}}(undef, model.nElements)

        for (i, element) in enumerate(model.elements)
            Ψ[i] = element.Ψ
            E[i] = element.section.E
            G[i] = element.section.G
            A[i] = element.section.A
            Ix[i] = element.section.Ix
            Iy[i] = element.section.Iy
            J[i] = element.section.J

            nodeids[i] = Asap.nodeids(element)
            dofids[i] = element.globalID
        end

        #Variable processing
        vals, lowerbounds, upperbounds = process_variables!(variables)

        #Indexer
        indexer = FrameOptIndexer(variables)

        #topology
        C = Asap.connectivity(model)
        freeids = model.freeDOFs

        #loads
        P = model.P
        Pf = model.Pf

        #sparsity pattern of K
        inzs = all_inz(model)
        cp = model.S.colptr
        rv = model.S.rowval
        nnz = length(model.S.nzval)

        #element releases
        releases = get_element_release.(model.elements)


        new{T}(
            model,
            vals,
            indexer,
            variables,
            releases,
            X,
            Y,
            Z,
            Ψ,
            E,
            G,
            A,
            Ix,
            Iy,
            J,
            P,
            Pf,
            C,
            lowerbounds,
            upperbounds,
            cp,
            rv,
            nnz,
            inzs,
            freeids,
            nodeids,
            dofids,
            model.nDOFs
        )


    end

end