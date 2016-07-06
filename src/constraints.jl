function populateConstraints(ordgdp::ORDGDP, p::problemData, scen_idx::Int64)

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
 
    
    
    # VOLTAGE CONSTRAINTS
    for i = 1:numEdges
      edge = p.EDGES[i]
      
        # TODO These references should be const
        idx1 = p.hashTableNodes[edge.node1id]
        idx2 = p.hashTableNodes[edge.node2id]
        lineCodeIdx = p.hashTableLineCodes[edge.lineCode] 
        rmatrix = p.LINECODES[lineCodeIdx].rmatrix * edge.length
        xmatrix = p.LINECODES[lineCodeIdx].xmatrix * edge.length
        
        if p.EDGES[i].numPhases == 1
            for j = 1:numPhases
                if p.EDGES[i].hasPhase[j]
                    @constraint(mip_model, voltageVariable[idx2,j] - voltageVariable[idx1,j] + voltageOffsetVariable[i,j] + 2.0*rmatrix[1,j,j]*flowRealVariable[i,j] + 2.0*xmatrix[1,j,j]*flowReactiveVariable[i,j]== 0.0)
                end
            end
        else
            for j = 1:numPhases 
                @constraint(mip_model, voltageVariable[idx2,j] - voltageVariable[idx1,j] + voltageOffsetVariable[i,j] + sum{2.0*rmatrix[mod(numPhases-j+l, numPhases)+1,j,l]*flowRealVariable[i,l] + 2.0*xmatrix[mod(numPhases-j+l, numPhases)+1,j,l]*flowReactiveVariable[i,l], l = 1:numPhases} == 0.0)              
            end
        end
    end
  
    
    # TODO Properly calculate this constant using p.u. or actual max voltage
    voltageMAX = 1e+5
    # VOLTAGE OFFSET CONSTRAINT
    for j = 1:numEdges
        for k = 1:numPhases
            if p.EDGES[j].hasPhase[k]
                @constraint(mip_model, voltageOffsetVariable[j,k] <= voltageMAX - lineExistsVariable[j] * voltageMAX)
                @constraint(mip_model, voltageOffsetVariable[j,k] >= -voltageMAX + lineExistsVariable[j] * voltageMAX)
            end
        end
    end

    # LINE EXISTS CONSTRAINT
    for j = 1:numEdges
        @constraint(mip_model, lineExistsVariable[j] - lineUseVariable[j] + switchUseVariable[j] == 0.0)
    end

    # VOLTAGE LIMIT CONSTRAINT
    for j = 1:numNodes
        refVoltage = p.NODES[j].refVoltage
        maxVoltage = p.NODES[j].maxVoltage
        minVoltage = p.NODES[j].minVoltage
        for k = 1:numPhases
            if p.NODES[j].hasPhase[k]
                @constraint(mip_model, voltageVariable[j,k] <= maxVoltage^2 * 1e+3)
                @constraint(mip_model, voltageVariable[j,k] >= minVoltage^2 * 1e+3)
            end
        end
    end

    # LINE CONSTRAINT
    for j = 1:numEdges
        for k = 1:numPhases
            if p.EDGES[j].hasPhase[k]
                @constraint(mip_model, flowRealVariable[j,k] <= p.EDGES[j].capacity * lineDirectionVariableForward[j])
                @constraint(mip_model, flowRealVariable[j,k] >= -p.EDGES[j].capacity * lineDirectionVariableBackward[j])

                @constraint(mip_model, flowReactiveVariable[j,k] <= p.EDGES[j].capacity * lineDirectionVariableForward[j])
                @constraint(mip_model, flowReactiveVariable[j,k] >= -p.EDGES[j].capacity * lineDirectionVariableBackward[j])
            end
        end
    end

    # LINE DIRECTION CONSTRAINT
    for j = 1:numEdges
        @constraint(mip_model, lineExistsVariable[j] == lineDirectionVariableForward[j] + lineDirectionVariableBackward[j]) 
    end 

    PhaseVariation = p.phaseVariation 
    # PHASE VARIATION CONSTRAINT
    for j = 1:numEdges
        if p.EDGES[j].numPhases > 1 && p.EDGES[j].isTransformer
            for k = 1:numPhases
                @constraint(mip_model, flowRealVariable[j,k] * (numPhases + PhaseVariation - 1.0)/numPhases + sum{flowRealVariable[j,l] * (PhaseVariation - 1.0)/numPhases, l = 1:numPhases} >= 0)
                @constraint(mip_model, flowRealVariable[j,k] * (PhaseVariation + 1.0 - numPhases)/numPhases + sum{flowRealVariable[j,l] * (PhaseVariation + 1.0)/numPhases, l = 1:numPhases} >= 0)

                @constraint(mip_model, flowReactiveVariable[j,k] * (numPhases + PhaseVariation - 1.0)/numPhases + sum{flowReactiveVariable[j,l] * (PhaseVariation - 1.0)/numPhases, l = 1:numPhases} >= 0)
                @constraint(mip_model, flowReactiveVariable[j,k] * (PhaseVariation + 1.0 - numPhases)/numPhases + sum{flowReactiveVariable[j,l] * (PhaseVariation + 1.0)/numPhases, l = 1:numPhases} >= 0)
            end
        end
    end

    # EVERYDAY OPERATION = NO DAMAGE
    # EMRE, i changed this conditional.  Is it what you meant?    
    #if scen_idx == 0  
    if scen_idx != 0
        # HARDEN CONSTRAINT        
        # seem right EMRE???
        for j = 1:numEdges
          edgeid = p.EDGES[j].id
          is_harden_disable = false
          if haskey(p.SCENARIOS[scen_idx].harden_disabled_edges,edgeid)
            is_hard_disable = p.SCENARIOS[scen_idx].harden_disabled_edges[edge_id]
          end
          
          if is_harden_disable
            @constraint(mip_model, lineHardenVariable[j] <= 0.0)
          else
            @constraint(mip_model, lineHardenVariable[j] <= 1.0)
          end
        end
#        for j = 1:numEdges
 #           if p.HARDENED_DISABLED[j].data
  #              @constraint(mip_model, lineHardenVariable[j] <= 0.0)
   #         else
    #            @constraint(mip_model, lineHardenVariable[j] <= 1.0)
     #       end
     #   end
        
        # DAMAGE CONSTRAINT
        supportable = zeros(Bool, numEdges)
        for j = 1:length(p.HARDEN_COST)
            idx = hashTableEdges[p.HARDEN_COST[j].id]
            supportable[idx] = true
        end
        
        # EMRE, seem right????
        for j = 1:numEdges
          edgeid = p.EDGES[j].id
          is_disable = false
          if haskey(p.SCENARIOS[scen_idx].disabled_edges,edgeid)
            is_disable = p.SCENARIOS[scen_idx].disabled_edges[edge_id]
          end
          
          if disable
            if supportable[j]
              @constraint(mip_model, lineUseVariable[j] - lineHardenVariable[j]== 0)
            else
              @constraint(mip_model, lineUseVariable[j] == 0)
            end
          end
        end        
        #for j = 1:length(p.DISABLED)
         #   idx = hashTableEdges[p.DISABLED[j].id]
          #  if p.DISABLED[j].data
           #     if supportable[idx]
            #        @constraint(mip_model, lineUseVariable[idx] - lineHardenVariable[idx]== 0)
             #   else
              #      @constraint(mip_model, lineUseVariable[idx] == 0)
              #  end
         #   end
        #end
        
        
    end

    # LOAD CONSTRAINT
    for j = 1:numNodes
        for k = 1:numPhases
            if p.NODES[j].hasPhase[k]
                @constraint(mip_model, loadRealVariable[j,k] == sum{loadServeVariable[l] * p.LOADS[l].maxRealPhase[k], l in p.NODES[j].LoadList})  
                @constraint(mip_model, loadReactiveVariable[j,k] == sum{loadServeVariable[l] * p.LOADS[l].maxReactivePhase[k], l in p.NODES[j].LoadList})  
            end
        end
    end

    # GENERATION CONSTRAINT
    maxRealPhase = zeros(Float64, numNodes, numPhases)
    maxReactivePhase = zeros(Float64, numNodes, numPhases)
    
    for j = 1:numNodes
        for l in p.NODES[j].GeneratorList
            for k = 1:numPhases
              maxRealPhase[j,k] += p.GENERATORS[l].maxRealPhase[k] 
              maxReactivePhase[j,k] += p.GENERATORS[l].maxReactivePhase[k] 
            end
        end
    end
    maxAdded = zeros(Float64, numGenerators)
    for j = 1:length(p.MAX_MICROGRID)
        idx = p.hashTableGenerators[p.MAX_MICROGRID[j].id]
        maxAdded[idx] = p.MAX_MICROGRID[j].data
    end
    for j = 1:numNodes
        for k = 1:numPhases
            @constraint(mip_model, generatorRealVariable[j,k] <= maxRealPhase[j,k] + sum{facilityVariable[l] * maxAdded[l], l in p.NODES[j].GeneratorList})
            @constraint(mip_model, generatorReactiveVariable[j,k] <= maxReactivePhase[j,k] + sum{facilityVariable[l] * maxAdded[l], l in p.NODES[j].GeneratorList})
        end
    end
    
    # FLOW BALANCE CONSTRAINT
    for j = 1:numNodes
        for k = 1:numPhases
            if p.NODES[j].hasPhase[k] 
                @constraint(mip_model, generatorRealVariable[j,k] - loadRealVariable[j,k] == sum{flowRealVariable[l,k], l in p.NODES[j].EdgeOutList} - sum{flowRealVariable[l,k], l in p.NODES[j].EdgeInList})             
                @constraint(mip_model, generatorReactiveVariable[j,k] - loadReactiveVariable[j,k] == sum{flowReactiveVariable[l,k], l in p.NODES[j].EdgeOutList} - sum{flowReactiveVariable[l,k], l in p.NODES[j].EdgeInList})             
            end
        end
    end

    # CYCLE ELIMINATION CONSTRAINT
    for j = 1:length(p.CYCLES)
        cycle = p.CYCLES[j]
        EdgeList = Int64[]
        C = length(cycle)
        for l = 1:numEdges
            idx1 = p.hashTableNodes[p.EDGES[l].node1id]
            idx2 = p.hashTableNodes[p.EDGES[l].node2id]
            if in(idx1, cycle) && in(idx2, cycle)
                idx = p.hashTableUniqueEdges[l]
                if !in(idx, EdgeList)
                    push!(EdgeList, idx)
                end
            end
        end
        @constraint(mip_model, sum{lineCycleVariable[l] - switchCycleVariable[l], l in EdgeList} <= C - 1.0)
    end
    for j = 1:numEdges
        @constraint(mip_model, lineUseVariable[j] <= lineCycleVariable[hashTableUniqueEdge[j]])
        @constraint(mip_model, switchCycleVariable[hashTableUniqueEdge[j]] <= lineCycleVariable[hashTableUniqueEdge[j]])
    end
    
end

function populateTyingConstraints(ordgdp::ORDGDP, p::problemData)

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
 
    subVariable = Any[]
    push!(subVariable, lineUseVariable)
    push!(subVariable, switchUseVariable)
    push!(subVariable, lineHardenVariable)
    push!(subVariable, facilityVariable)
    

    for i = 1:length(ordgdp.masterVariable)
        for j = 1:length(ordgdp.masterVariable[i])
            @constraint(mip_model, ordgdp.masterVariable[i][j] >= subVariable[i][j])
        end
    end

end
