"""
    NetworkOptParams(model::Network, variables::Vector{NetworkVariable})

Contains all information and fields necessary for optimization.
"""
struct NetworkOptParams <: AbstractOptParams
    network::Network #the reference truss model for optimization
    values::Vector{Float64} #design variables
    indexer::NetworkOptIndexer #pointers to design variables and full variables
    variables::Vector{AbstractVariable}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    q::Vector{Float64} #all element young's modulii |n_element|
    Pn::Matrix{Float64}
    C::SparseMatrixCSC{Int64, Int64} #connectivity matrix
    Cn::SparseMatrixCSC{Int64, Int64}
    Cf::SparseMatrixCSC{Int64, Int64}
    lb::Vector{Float64} #lower bounds of variables
    ub::Vector{Float64} #upper bounds of variables
    N::Vector{Int64}
    F::Vector{Int64}

    function NetworkOptParams(network::Network, variables::Vector{T}) where T <: NetworkVariable
        @assert network.processed "network must be processed"

        #extract global parameters
        xyz = network.xyz
        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]
        q = network.q

        #Variable processing
        vals, lowerbounds, upperbounds = process_variables!(variables)

        #Indexer
        indexer = NetworkOptIndexer(variables)

        #generate a truss optimization problem
        new(network,
            vals,
            indexer,
            variables,
            X,
            Y,
            Z,
            q,
            network.Pn,
            network.C,
            network.Cn,
            network.Cf,
            lowerbounds,
            upperbounds,
            network.N,
            network.F
            )

    end
end