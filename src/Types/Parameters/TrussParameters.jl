"""
    TrussOptParams(model::TrussModel, variables::Vector{TrussVariable})

Contains all information and fields necessary for optimization.
"""
struct TrussOptParams <: AbstractOptParams
    model::TrussModel #the reference truss model for optimization
    values::Vector{Float64} #design variables
    indexer::TrussOptIndexer #pointers to design variables and full variables
    variables::Vector{AbstractVariable}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    E::Vector{Float64} #all element young's modulii |n_element|
    A::Vector{Float64} #all element areas |n_element|
    P::Vector{Float64} # External load vector
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

    function TrussOptParams(model::TrussModel, variables::Vector{T}) where T <: TrussVariable
        
        #model must be pre-proces
        model.processed || (Asap.process!(model))

        #extract global parameters
        xyz = node_positions(model)
        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]
        E = getproperty.(getproperty.(model.elements, :section), :E)
        A = getproperty.(getproperty.(model.elements, :section), :A)

        #Variable processing
        vals, lowerbounds, upperbounds = process_variables!(variables)

        #Indexer
        indexer = TrussOptIndexer(variables)

        #topology
        nodeids = Asap.nodeids.(model.elements)
        dofids = getproperty.(model.elements, :globalID)
        C = Asap.connectivity(model)
        freeids = model.freeDOFs

        #external load
        P = model.P

        #sparsity pattern of K
        inzs = all_inz(model)
        cp = model.S.colptr
        rv = model.S.rowval
        nnz = length(model.S.nzval)

        #generate a truss optimization problem
        new(model, 
            vals, 
            indexer, 
            variables, 
            X, 
            Y, 
            Z, 
            E, 
            A, 
            P,
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