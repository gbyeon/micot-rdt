

function populateVariables(ordgdp::ORDGDP, p::problemData, scen_idx::Int64, recordMaster::Bool)

    numPhases = p.numPhases
    numNodes = length(p.NODES)
    numEdges = length(p.EDGES)
    numGenerators = length(p.GENERATORS)
    numLoads = length(p.LOADS)

    # GET MODEL REFERENCE
    mip_model = ordgdp.mip_model
 
    # EDGE BASED VARIABLES
    ordgdp.lineCycleVariable = @variable(mip_model, lineCycleVariable[1:numEdges], Bin)
    ordgdp.switchCycleVariable = @variable(mip_model, switchCycleVariable[1:numEdges], Bin)
    ordgdp.lineUseVariable = @variable(mip_model, lineUseVariable[1:numEdges], Bin)
    ordgdp.lineExistsVariable = @variable(mip_model, lineExistsVariable[1:numEdges], Bin)
    ordgdp.lineDirectionVariableForward = @variable(mip_model, lineDirectionVariableForward[1:numEdges], Bin)
    ordgdp.lineDirectionVariableBackward = @variable(mip_model, lineDirectionVariableBackward[1:numEdges], Bin)
    ordgdp.switchUseVariable = @variable(mip_model, switchUseVariable[1:numEdges], Bin)
    ordgdp.lineHardenVariable = @variable(mip_model, lineHardenVariable[1:numEdges], Bin)
    ordgdp.flowRealVariable = @variable(mip_model, flowRealVariableA[1:numEdges,1:numPhases])
    ordgdp.flowReactiveVariable = @variable(mip_model, flowReactiveVariableA[1:numEdges,1:numPhases])
    ordgdp.voltageOffsetVariable = @variable(mip_model, voltageOffsetVariableA[1:numEdges,1:numPhases])

    # NODE BASED VARIABLES
    ordgdp.voltageVariable = @variable(mip_model, voltageVariableA[1:numNodes,1:numPhases])
    ordgdp.loadRealVariable = @variable(mip_model, loadRealVariableA[1:numNodes,1:numPhases])
    ordgdp.loadReactiveVariable = @variable(mip_model, loadReactiveVariableA[1:numNodes,1:numPhases])
    ordgdp.generatorRealVariable = @variable(mip_model, generatorRealVariableA[1:numNodes,1:numPhases])
    ordgdp.generatorReactiveVariable = @variable(mip_model, generatorReactiveVariableA[1:numNodes,1:numPhases])

    # LOAD BASED VARIABLES
    ordgdp.loadServeVariable = @variable(mip_model, loadServeVariable[1:numLoads], Bin)

    # GENERATOR BASED VARIABLES
    ordgdp.facilityVariable = @variable(mip_model, facilityVariable[1:numGenerators], Bin)

    # RECORD MASTER VARIABLES
    # TODO make it user defined!!!
    if recordMaster
        masterVariable = Any[]
        push!(masterVariable, lineUseVariable)
        push!(masterVariable, switchUseVariable)    
        push!(masterVariable, lineHardenVariable)    
        push!(masterVariable, facilityVariable)   
        ordgdp.masterVariable = masterVariable
    end

end
