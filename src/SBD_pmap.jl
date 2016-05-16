# Following code performs a parallel Scenario-based Decomposition (SBD)
# of an arbitrary two stage mixed-integer program where only first
# stage variables appear at the objective function.
#
# User must provide a constructor function such that:
# function constructor(master_model, master_var, max_ind)
#   if master_model == []
#   // CREATE SUBPROBLEM max_ind
#   // ADD CONSTRAINTS S.T. sub_var[i][j] == master_var[i][j]
#   // SOLVE SUBPROBLEM
#   // RETURN OBJECTIVEVALUE, SOLUTION CORRESPONDING TO sub_var
#   else
#   // APPEND SUBPROBLEM max_ind TO master_model
#   end
# end
#
# if user wants to solve the master problem with an exact method, please wrap
# solve from JuMP with a new method and provide it as the master_solver

using JuMP

function getObjective(cost, curr_solution)

    currObj = 0.0
    for j = 1:length(curr_solution)
        for k = 1:length(curr_solution[j])
            currObj += cost[j][k] * curr_solution[j][k]
        end
    end

    return currObj
end

# Takes the max of the curr_solution, over scenarios
function Greedy(curr_solution)

    num_sub_prob = length(curr_solution)
    greedy_solution = Any[]
    for j = 1:length(curr_solution[1])
        push!(greedy_solution, zeros(length(curr_solution[1][j])))
    end

    for j = 1:length(curr_solution[1])   
        for k = 1:length(curr_solution[1][j])
            max_val = 0.0
            for i = 1:num_sub_prob
                if max_val < curr_solution[i][j][k]
                    max_val = curr_solution[i][j][k]
                end
            end
            greedy_solution[j][k] = max_val
        end
    end
    return greedy_solution

end

function findMax(arr)
    max_val = -Inf
    max_ind = -1
    for i = 1:length(arr)
        if arr[i] > max_val
            max_val = arr[i]
            max_ind = i
        end
    end
    return max_val, max_ind
end

function solveSubProblems(constructor::Function, master_solution, curr_objective, curr_solution, num_sub_prob, fixVals::Bool)

    output = pmap((a1,a2,a3,a4)->constructor(a1,a2,a3,a4), [[] for i = 1:num_sub_prob], [fixVals for i = 1:num_sub_prob], [master_solution for i = 1:num_sub_prob], 1:num_sub_prob)
    @show output
    # TODO collect only subproblem solutions that were not considered
    # record subproblem solutions
    for i = 1:num_sub_prob
        curr_objective[i] = output[i][1]
        for j = 1:length(curr_solution[i])
            for k = 1:length(curr_solution[i][j])
                curr_solution[i][j][k] = output[i][2][j][k]
            end
        end
    end

    # perform greedy
    greedy_solution = Greedy(curr_solution)

    return greedy_solution

end

# master_var is array of array of array
function ScenarioBasedDecomposition_pmap(constructor::Function, master_solver::Function, cost, master_object, master_var, num_sub_prob)

    println("SBD STARTED")
 
    start = time()
   
    OPTIMALITY_TOL = 1e-4
    z_LB = 0
    addedScenarios = zeros(Int, num_sub_prob)
    
    num_master_var = length(master_var)

    curr_objective = zeros(num_sub_prob)
    # create lagrange array 
    curr_solution = Any[]
    for i = 1:num_sub_prob
        curr_sub = Any[]
        for j = 1:num_master_var
            push!(curr_sub, zeros(length(master_var[j])))
        end
        push!(curr_solution, curr_sub)
    end

    master_solution = Any[]
    for j = 1:num_master_var
        push!(master_solution, zeros(length(master_var[j])))
    end

    # perform initial greedy
    greedy_solution_ini = solveSubProblems(constructor, master_solution, curr_objective, curr_solution, num_sub_prob, false)
    greedy_solution = copy(greedy_solution_ini)
    z_UB = getObjective(cost, greedy_solution)

    @show greedy_solution_ini 
    # MAIN LOOP
    # iterate until greedy objective is lower than master objective
    # with a given tolerance
    iter = 0
    while (z_UB - z_LB) > OPTIMALITY_TOL

        # set initial descent from best known solution
        for i = 1:num_master_var
            for j = 1:length(master_var[i])
                master_solution[i][j] = greedy_solution[i][j]
            end
        end

        # solve master problem
        master_solver(master_object.mip_model, master_solution, master_var)
        z_LB = getObjective(cost, master_solution)

        # TODO activate the following portion after master solve
        for k = 1:num_sub_prob
            if addedScenarios[k] == 1
                for i = 1:num_master_var
                    for j = 1:length(master_var[i])
                        curr_solution[k][i][j] = master_solution[i][j]
                    end
                end
            else
                for i = 1:num_master_var
                    for j = 1:length(master_var[i])
                        curr_solution[k][i][j] = greedy_solution_ini[i][j]
                    end
                end
            end
        end

        greedy_solution = Greedy(curr_solution)

        solveSubProblems(constructor, master_solution, curr_objective, curr_solution, num_sub_prob, true)
        z_UB = getObjective(cost, greedy_solution)

        @show greedy_solution
        @show z_LB, z_UB
        # find max violated objective
        max_val, max_ind = findMax(curr_objective)
        if abs(max_val) > OPTIMALITY_TOL
            constructor(master_object, false, master_var, max_ind)
            addedScenarios[max_ind] = 1
            iter += 1
        else 
            # solution optimum
            break
        end
    end
    return z_LB
end
