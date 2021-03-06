function populateInfeasibilityMinimization(ordgdp::ORDGDP)

    numPhases = p.numPhases
    numNodes = length(p.NODES)
    numEdges = length(p.EDGES)
    numGenerators = length(p.GENERATORS)
    numLoads = length(p.LOADS)

    # GET REFERENCES
    mip_model = ordgdp.mip_model

    lineCycleVariable = ordgdp.lineCycleVariable
    switchCycleVariable = ordgdp.switchCycleVariable
    lineUseVariable = ordgdp.lineUseVariable
    lineExistsVariable = ordgdp.lineExistsVariable
    lineDirectionVariableForward = ordgdp.lineDirectionVariableForward
    lineDirectionVariableBackward = ordgdp.lineDirectionVariableBackward
    switchUseVariable = ordgdp.switchUseVariable
    lineHardenVariable = ordgdp.lineHardenVariable
    flowRealVariable = ordgdp.flowRealVariable
    flowReactiveVariable = ordgdp.flowReactiveVariable
    voltageOffsetVariable = ordgdp.voltageOffsetVariable
    voltageVariable = ordgdp.voltageVariable
    loadRealVariable = ordgdp.loadRealVariable
    loadReactiveVariable = ordgdp.loadReactiveVariable
    generatorRealVariable = ordgdp.generatorRealVariable
    generatorReactiveVariable = ordgdp.generatorReactiveVariable
    loadServeVariable = ordgdp.loadServeVariable
    facilityVariable = ordgdp.facilityVariable
 
    @variable(mip_model, infCriticalRealVariable[1:numPhases] >= 0)
    @variable(mip_model, infCriticalReactiveVariable[1:numPhases] >= 0)

    @variable(mip_model, infTotalRealVariable[1:numPhases] >= 0)
    @variable(mip_model, infTotalReactiveVariable[1:numPhases] >= 0)

    CriticalLoadList = Int64[]
    for j = 1:length(p.IS_CRITICAL_LOAD)
        idx = p.hashTableLoads[p.IS_CRITICAL_LOAD[j].id]
        if p.IS_CRITICAL_LOAD[j]
            push!(CriticalLoadList, idx)
        end
    end
    for k = 1:numPhases
        @constraint(mip_model, sum{loadRealVariable[j,k], j in CriticalLoadList} + infCriticalRealVariable[k] >= p.CriticalRealDemand[k])
        @constraint(mip_model, sum{loadReactiveVariable[j,k], j in CriticalLoadList} + infCriticalReactiveVariable[k] >= p.CriticalReactiveDemand[k])

        @constraint(mip_model, sum{loadRealVariable[j,k], j in 1:numLoads} + infTotalRealVariable[k] >= p.TotalRealDemand[k])
        @constraint(mip_model, sum{loadReactiveVariable[j,k], j in 1:numLoads} + infTotalReactiveVariable[k] >= p.TotalReactiveDemand[k])
    end

    @objective(mip_model, Min, sum{infCriticalRealVariable[k] + infCriticalReactiveVariable[k] + infTotalRealVariable[k] + infTotalReactiveVariable[k], k in 1:numPhases})

end

function populateMinimization(ordgdp::ORDGDP, p::problemData)

    numPhases = p.numPhases
    numNodes = length(p.NODES)
    numEdges = length(p.EDGES)
    numGenerators = length(p.GENERATORS)
    numLoads = length(p.LOADS)

    # GET REFERENCES
    mip_model = ordgdp.mip_model

    lineCycleVariable = ordgdp.lineCycleVariable
    switchCycleVariable = ordgdp.switchCycleVariable
    lineUseVariable = ordgdp.lineUseVariable
    lineExistsVariable = ordgdp.lineExistsVariable
    lineDirectionVariableForward = ordgdp.lineDirectionVariableForward
    lineDirectionVariableBackward = ordgdp.lineDirectionVariableBackward
    switchUseVariable = ordgdp.switchUseVariable
    lineHardenVariable = ordgdp.lineHardenVariable
    flowRealVariable = ordgdp.flowRealVariable
    flowReactiveVariable = ordgdp.flowReactiveVariable
    voltageOffsetVariable = ordgdp.voltageOffsetVariable
    voltageVariable = ordgdp.voltageVariable
    loadRealVariable = ordgdp.loadRealVariable
    loadReactiveVariable = ordgdp.loadReactiveVariable
    generatorRealVariable = ordgdp.generatorRealVariable
    generatorReactiveVariable = ordgdp.generatorReactiveVariable
    loadServeVariable = ordgdp.loadServeVariable
    facilityVariable = ordgdp.facilityVariable
 
    CriticalLoadList = Int64[]
    for j = 1:length(p.IS_CRITICAL_LOAD)
        idx = p.hashTableLoads[p.IS_CRITICAL_LOAD[j].id]
        if p.IS_CRITICAL_LOAD[j].data
            push!(CriticalLoadList, idx)
        end
    end
    
    
    
    # TODO RBENT - Another place for refactoring... constraints are living here
    @constraint(mip_model, sum{loadRealVariable[j,k], j in CriticalLoadList, k in 1:numPhases} >= p.CriticalRealDemand)
    @constraint(mip_model, sum{loadReactiveVariable[j,k], j in CriticalLoadList,k in 1:numPhases} >= p.CriticalReactiveDemand)
    @constraint(mip_model, sum{loadRealVariable[j,k], j in 1:numLoads, k in 1:numPhases} >= p.TotalRealDemand)
    @constraint(mip_model, sum{loadReactiveVariable[j,k], j in 1:numLoads, k in 1:numPhases} >= p.TotalReactiveDemand)
           
#    for k = 1:numPhases
 #       @constraint(mip_model, sum{loadRealVariable[j,k], j in CriticalLoadList} >= p.CriticalRealDemand[k])
  #      @constraint(mip_model, sum{loadReactiveVariable[j,k], j in CriticalLoadList} >= p.CriticalReactiveDemand[k])
   #     @constraint(mip_model, sum{loadRealVariable[j,k], j in 1:numLoads} >= p.TotalRealDemand[k])
    #    @constraint(mip_model, sum{loadReactiveVariable[j,k], j in 1:numLoads} >= p.TotalReactiveDemand[k])
    #end

    maxAdded = zeros(Float64, numGenerators)
    for j = 1:length(p.MAX_MICROGRID)
        idx = p.hashTableGenerators[p.MAX_MICROGRID[j].id]
        maxAdded[idx] = p.MAX_MICROGRID[j].data
    end

    FACILITY_COST = zeros(Float64, numGenerators)
    for j = 1:length(p.MICROGRID_COST)
        idx = p.hashTableGenerators[p.MICROGRID_COST[j].id]
        FACILITY_COST[idx] += maxAdded[idx] * p.MICROGRID_COST[j].data
    end
    for j = 1:length(p.MICROGRID_FIXED_COST)
        idx = p.hashTableGenerators[p.MICROGRID_FIXED_COST[j].id]
        FACILITY_COST[idx] += p.MICROGRID_FIXED_COST[j].data
    end
    
    #TODO, potential bug here, we might be able to turn off the line without a switch here here.
    #TODO, RBENT need to only post variables for switch cost and line construction cost
    
    # JUMP doesn't like Inf in coefficients, so just doing a big number until we fix the above thing
    LINE_SWITCH_COST = zeros(Float64,length(p.LINE_SWITCH_COST))
    for j=1:length(p.LINE_SWITCH_COST)
      if p.LINE_SWITCH_COST[j].data == Inf
        LINE_SWITCH_COST[j] = 1e10
      else
        LINE_SWITCH_COST[j] = p.LINE_SWITCH_COST[j].data
      end
    end
    
    HARDEN_COST = zeros(Float64,length(p.HARDEN_COST))
    for j=1:length(p.HARDEN_COST)
      if p.HARDEN_COST[j].data == Inf
        HARDEN_COST[j] = 1e10
      else
        HARDEN_COST[j] = p.HARDEN_COST[j].data
      end
    end
    
    
    @objective(mip_model, Min, 
        sum{p.LINE_CONSTRUCTION_COST[j].data * lineUseVariable[p.hashTableEdges[p.LINE_CONSTRUCTION_COST[j].id]], j in 1:length(p.LINE_CONSTRUCTION_COST)} +
        sum{LINE_SWITCH_COST[j] * switchUseVariable[p.hashTableEdges[p.LINE_SWITCH_COST[j].id]], j in 1:length(p.LINE_SWITCH_COST)} +
        sum{HARDEN_COST[j] * lineHardenVariable[p.hashTableEdges[p.HARDEN_COST[j].id]], j in 1:length(p.HARDEN_COST)} +
        sum{FACILITY_COST[j] * facilityVariable[j], j in 1:numGenerators})
        
end

