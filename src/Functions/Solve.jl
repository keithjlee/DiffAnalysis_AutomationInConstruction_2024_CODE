"""
    solve_u(K::SparseMatrixCSC{Float64, Int64}, p::TrussOptParams)

Displacement of free DOFs
"""
function solve_u(K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams, alg::LinearSolve.SciMLLinearSolveAlgorithm)
    linearproblem = LinearProblem(K[p.freeids, p.freeids], p.P[p.freeids])
    sol = LinearSolve.solve(linearproblem, alg)

    return sol.u
end

function solve_u_noadjoint(K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams, alg::LinearSolve.SciMLLinearSolveAlgorithm)
    linearproblem = LinearProblem(K[p.freeids, p.freeids], p.P[p.freeids])

    sol = LinearSolve.solve(linearproblem, alg)
    return sol.u
end

function ChainRulesCore.rrule(::typeof(solve_u), K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams, alg::LinearSolve.SciMLLinearSolveAlgorithm)
    linearproblem = LinearProblem(K[p.freeids, p.freeids], p.P[p.freeids])
    linsolve = LinearSolve.init(linearproblem, alg)
    sol1 = LinearSolve.solve!(linsolve)
    u = sol1.u

    function solve_u_pullback(ū)

        #initialize
        dudK = zeros(p.n, p.n)

        #sensitivities w/r/t active DOFs
        ku = linsolve.cacheval \ ū

        #assign to proper indices
        @views dudK[p.freeids, p.freeids] = ku * u'

        return NoTangent(), -dudK, NoTangent(), NoTangent()
    end

    return u, solve_u_pullback
end