using FactCheck

include("SBD_pmap.jl")
include("VNS.jl")

@everywhere using JuMP
@everywhere using CPLEX
@everywhere begin

type masterData
    mip_model
    function masterData(mip_model)
        m = new()
        m.mip_model = mip_model
        return m
    end
end

function fixSolution(x, values)

   for i in 1:length(x)
      setLower(x[i], round(values[i]))
      setUpper(x[i], round(values[i]))
   end

end

function generate3Knapsack(master_object, fixVals, master_var, scen_id)
    println("ENTERING SCENARIO FUNCTION")

    start = time()
 
    c = [2,3,1]
    if scen_id == 1
        a = [1,1,1]
        b = 1
    elseif scen_id == 2
        a = [3,1,1]
        b = 3
    elseif scen_id == 3
        a = [3,3,1]
        b = 6
    end
    println("DATA CREATED")
    if master_object != []
        curr_model = master_object.mip_model
    else
        curr_model = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0))
    end

    curr_solution = Any[]
    sub_solution = Any[]

    println("POPULATING SCENARIO $scen_id")
    @variable(curr_model, x[1:3], Bin)
    if master_object == []
        println("STARTED SCENARIO $scen_id")
        if fixVals
            fixSolution(x, master_var[1])
            @variable(curr_model, inf_var >= 0)
            @constraint(curr_model, vecdot(a, x) + inf_var >= b)
            @setObjective(curr_model, Min, inf_var)
        else
            @constraint(curr_model, vecdot(a, x) >= b)
            @setObjective(curr_model, Min, vecdot(c, x))
        end
        println("SOLVING SCENARIO $scen_id")
        solve(curr_model)
        println("SOLVED SCENARIO $scen_id")
        for i = 1:3
            push!(sub_solution, getValue(x[i]))
        end
        push!(curr_solution, sub_solution)

        println("process id $(myid()) took $(time()-start) seconds.")
        return getObjectiveValue(curr_model), curr_solution
    else
        @constraint(curr_model, vecdot(a, x) >= b)
        for j = 1:3
            @constraint(curr_model, master_var[1][j] >= x[j])
        end
    end
end

end

function generateMaster()
    master_var = Any[]
    master_model = Model(solver=CplexSolver()) 
    c = [2,3,1]
    @variable(master_model, x[1:3], Bin)
    @setObjective(master_model, Min, vecdot(c, x))
    push!(master_var, x)
    return master_model, master_var 
end

function master_solver(master_model, master_solution, master_var)
    lp_relaxation = Any[]
    push!(lp_relaxation, [0.5,0.5,0.5])
    return VariableNeighborhoodSearch(master_model, master_solution, lp_relaxation, master_var, Any[], Any[], 3600.0)
end

facts("Scenario-based Decomposition Test") do

    cost = Any[]
    c = [2,3,1]
    push!(cost, c)

    master_model, master_var = generateMaster()
    master_object = masterData(master_model)

    @fact ScenarioBasedDecomposition_pmap(generate3Knapsack, master_solver, cost, master_object, master_var, 3) --> 5.0
end
