using FactCheck
using JSON

@everywhere using JuMP
@everywhere using CPLEX
@everywhere begin

include("SBD_pmap.jl")
include("VNS.jl")
include("classes.jl")
include("objective.jl")
include("variables.jl")
include("constraints.jl")

mip_solver = CplexSolver()

# // -- GENERATE SCENARIO scen_idx ----------------------------
# TODO Move variables under their associated classes
function loadModel(ordgdp::ORDGDP, p::problemData, scen_idx::Int64, recordMaster::Bool)
   
    populateVariables(ordgdp, p, scen_idx, recordMaster)
    populateConstraints(ordgdp, p, scen_idx, recordMaster)

end

# // -- PARSER ------------------------------------------------

function loadProblemData(p::problemData)

end

function fixSolution(ordgdp, values)

    for i in 1:length(ordgdp.masterVariables)
        for j in 1:length(ordgdp.masterVariables[i])
            setLower(ordgdp.masterVariables[i][j], round(values[i][j]))
            setUpper(ordgdp.masterVariables[i][j], round(values[i][j]))
        end
    end

end

function generateScenario(master_ordgdp, fixVals, master_var, scen_idx)

    start = time()
 
    if master_model != []
        curr_ordgdp = master_ordgdp
    else
        problem_data = problemData()
        curr_ordgdp = ORDGDP(mip_solver, problem_data, scen_idx)
    end

    curr_solution = Any[]
    sub_solution = Any[]

    if master_model == []
        if fixVals
            fixSolution(curr_model, master_var)
            populateInfeasibilityMinimization(curr_model)
        else
            populateMinimization(curr_model)
        end
        solve(curr_model)
        for i = 1:length(curr_model.masterVariable)
            for j = 1:length(curr_model.masterVariable[i])
                push!(sub_solution, getValue(curr_model.masterVariable[i][j]))
            end
            push!(curr_solution, sub_solution)
        end

        return getObjectiveValue(curr_model), curr_solution
    else
        populateVariables(ordgdp, p, scen_idx, false)
        populateConstraints(ordgdp, p, scen_idx, false)
        populateTyingConstraints(ordgdp, p)
    end
end

end

function master_solver(master_model, master_solution, master_var)

    # MUST SOLVE LP RELAXATION FIRST
    lp_relaxation = Any[]

    return VariableNeighborhoodSearch(master_model, master_solution, lp_relaxation, master_var, Any[], Any[], 3600.0)
end

function solveORDGDP()

    problem_data = problemData()
    ordgdp = ORDGDP(mip_solver, problem_data, 0)

    master_ordgdp = generateScenario([], false, Any[], 0)

    #ScenarioBasedDecomposition_pmap(generateScenario, master_solver, cost, master_model, master_var, numSubProb)
end
