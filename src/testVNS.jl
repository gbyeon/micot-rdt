using JuMP
using FactCheck
using CPLEX

include("VNS.jl")

facts("Variable Neighborhood Search Test") do
    m = Model(solver=CplexSolver())

    @defVar(m, x >= 0, start = 1, Int)
    @defVar(m, z >= 0, start = 1, Int)
    @defVar(m, y >= 0, start = 1)

    @setObjective(m, Min, -3x - y + z)

    @addConstraint(m, 3x + 2y + z <= 20)

    @addConstraint(m, x + 3z <= 10)

    @addConstraint(m, z >= 1.5)

    master_var = Any[]
    push!(master_var, [x])
    push!(master_var, [z])

    incumbent = Any[]
    push!(incumbent, [1.0])
    push!(incumbent, [3.0])

    lp_relaxation = Any[]
    push!(lp_relaxation, [5.5])
    push!(lp_relaxation, [1.5])

    VariableNeighborhoodSearch(m, incumbent, lp_relaxation, master_var, Any[], Any[], 36000.0)

    @fact incumbent[1][1] --> roughly(4.0) 
    @fact incumbent[2][1] --> roughly(2.0) 

end
