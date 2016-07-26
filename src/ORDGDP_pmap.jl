using FactCheck
using JSON

@everywhere using JuMP
#@everywhere using CPLEX
#@everywhere using Cbc

@everywhere begin

include("SBD_pmap.jl")
include("VNS.jl")
include("classes.jl")
include("objective.jl")
include("variables.jl")
include("constraints.jl")

# Why is this a global variable.  TODO make this an input parameter
#mip_solver = CplexSolver()
#mip_solver = CbcSolver()

# // -- GENERATE SCENARIO scen_idx ----------------------------
# TODO Move variables under their associated classes
function loadModel(ordgdp::ORDGDP, p::problemData, scen_idx::Int64, recordMaster::Bool)
   
    populateVariables(ordgdp, p, scen_idx, recordMaster)
    populateConstraints(ordgdp, p, scen_idx)

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

# This is a kind of a multi-purposed generate scenario function
# If the master problem is empty, we solve the sub problem with master problem
# variables fixed. Otherwise, we just append this problem to the master problem
# This is designed to fit into the scenario-based decomposition constructor function
#
# master_ordgdp = the master problem, if empty we just solve the sub problem with fixVals
# fixed
#
# fixVals = the fixed variables of the master problem
# master_var = the master variables
# scen_idx = the scenario for which we are solving
# mip_solver = the solver for a mized integer problem
# TODO - RBENT - I think this should really be split into two seperate functions, they are doing 
# two very different things
function generateScenario(master_ordgdp, fixVals, master_var, scen_idx, mip_solver, problem_data)

    start = time()
 
    if master_ordgdp != []
        curr_ordgdp = master_ordgdp
    else
        #problem_data = problemData() RBENT -- I have added this as an input data
        curr_ordgdp = ORDGDP(mip_solver, problem_data, scen_idx, false)        
    end

    curr_solution = Any[]
    sub_solution = Any[]

    if master_ordgdp == []
        if fixVals
            fixSolution(curr_ordgdp, master_var)
            populateInfeasibilityMinimization(curr_ordgdp)
        else
            populateMinimization(curr_ordgdp, problem_data)
        end
        
        solve(curr_ordgdp.mip_model)
        
        for i = 1:length(curr_ordgdp.masterVariable)
            for j = 1:length(curr_ordgdp.masterVariable[i])
                push!(sub_solution, getValue(curr_model.masterVariable[i][j]))
            end
            push!(curr_solution, sub_solution)
        end

        return getobjectivevalue(curr_ordgdp.mip_model), curr_solution
    else
        populateVariables(ordgdp, p, scen_idx, false)
        populateConstraints(ordgdp, p, scen_idx, false)
        populateTyingConstraints(ordgdp, p)
    end
end

end

function master_solver(master_model, master_solution, master_var)

    # TODO - MUST SOLVE LP RELAXATION FIRST
    lp_relaxation = Any[]

    return VariableNeighborhoodSearch(master_model, master_solution, lp_relaxation, master_var, Any[], Any[], 3600.0)
end

function solveORDGDP(filename::AbstractString, mip_solver)

    # Read in all the data that we need to solve the problem
    problem_data = problemData(filename)
    
    # This creates the variables and constraints associated with the
    # base case
    ordgdp = ORDGDP(mip_solver, problem_data, 0, true)
        
    # Testing to see if the master problem solves, with no fixed variables
    master_ordgdp, sol = generateScenario([], false, Any[], 0, mip_solver, problem_data)

    println(sol)
    println(master_ordgdp) 
      
    #ScenarioBasedDecomposition_pmap(generateScenario, master_solver, cost, master_model, master_var, numSubProb)
end
