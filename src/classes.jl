# // -- ORDGDP CLASSES ---------------------------------

type dataNode
    id::AbstractString
    data
    function dataNode(name::AbstractString, x)
        d = new()
        d.id = name
        d.data = x
        return d
    end
end

type nodeData
    id::AbstractString
    x::Float64
    y::Float64
    EdgeInList::Vector{Int64}
    EdgeOutList::Vector{Int64}
    GeneratorList::Vector{Int64}
    LoadList::Vector{Int64}
    minVoltage::Float64
    maxVoltage::Float64
    refVoltage::Float64
    hasPhase::Vector{Bool}
    demand::Vector{Float64}
    hasGenerator::Bool
    function nodeData(id, x, y, minVoltage, maxVoltage, refVoltage, hasPhase, demand, hasGenerator)
        n = new()
        n.id = id
        n.x = x
        n.y = y
        n.EdgeInList = Int64[]
        n.EdgeOutList = Int64[]
        n.GeneratorList = Int64[]
        n.LoadList = Int64[]
        n.minVoltage = minVoltage
        n.maxVoltage = maxVoltage
        n.refVoltage = refVoltage
        n.hasPhase = hasPhase
        n.demand = demand
        n.hasGenerator = hasGenerator
        return n
    end
end

type generatorData
    id::AbstractString
    node_id::AbstractString
    hasPhase::Vector{Bool}
    maxRealPhase::Vector{Float64}
    maxReactivePhase::Vector{Float64}
    function generatorData(id, node_id, hasPhase, maxRealPhase, maxReactivePhase)
        g = new()
        g.id = id
        g.node_id = node_id
        g.hasPhase = hasPhase
        g.maxRealPhase = maxRealPhase
        g.maxReactivePhase = maxReactivePhase
        return g
    end
end

type loadData
    id::AbstractString
    node_id::AbstractString
    hasPhase::Vector{Bool}
    maxRealPhase::Vector{Float64}
    maxReactivePhase::Vector{Float64}
    function loadData(id, node_id, hasPhase, maxRealPhase, maxReactivePhase)
        l = new()
        l.id = id
        l.node_id = node_id
        l.hasPhase = hasPhase
        l.maxRealPhase = maxRealPhase
        l.maxReactivePhase = maxReactivePhase
        return l
    end   
end

type edgeData
    id::AbstractString
    node1id::AbstractString
    node2id::AbstractString
    hasPhase::Vector{Bool}
    capacity::Float64
    length::Float64
    numPhases::Float64
    isTransformer::Bool
    lineCode::Int64
    numPoles::Int64
    function edgeData(id, node1id, node2id, hasPhase, capacity, length, numPhases, isTransformer, lineCode, numPoles)
        e = new()
        e.id = id
        e.node1id = node1id
        e.node2id = node2id
        e.hasPhase = hasPhase
        e.capacity = capacity
        e.length = length
        e.numPhases = numPhases
        e.isTransformer = isTransformer
        e.lineCode = lineCode
        e.numPoles = numPoles
        return e
    end
end

type lineCodeData
    lineCode::Int64
    numPhases::Int64
    rmatrix
    xmatrix
    function lineCodeData(lineCode, numPhases, rmatrix, xmatrix)
        l = new()
        l.lineCode = lineCode
        l.numPhases = numPhases
        l.rmatrix = rmatrix
        l.xmatrix = xmatrix
        return l
    end 
end

type problemData
    # DATA
    # Power flow data
    numPhases::Int64    

    # Graph data
    NODES::Vector{nodeData}
    EDGES::Vector{edgeData}
    GENERATORS::Vector{generatorData}
    LOADS::Vector{loadData}
    LINECODES::Vector{lineCodeData}
    CYCLES

    # Hash tables
    hashTableNodes::Dict{AbstractString,Int64}
    hashTableEdges::Dict{AbstractString,Int64}
    hashTableUniqueEdges::Dict{AbstractString,Int64}
    hashTableGenerators::Dict{AbstractString,Int64}
    hashTableLoads::Dict{AbstractString,Int64}
    hashTableLineCodes::Dict{AbstractString,Int64}
    
    # Upgrade data
    LINE_CONSTRUCTION_COST::Vector{dataNode}
    MICROGRID_COST::Vector{dataNode}
    MICROGRID_FIXED_COST::Vector{dataNode}
    MAX_MICROGRID::Vector{dataNode}
    IS_CRITICAL_LOAD::Vector{dataNode}
    HARDEN_COST::Vector{dataNode}
    LINE_SWITCH_COST::Vector{dataNode}

    # Demand data
    CriticalRealDemand::Float64
    CriticalReactiveDemand::Float64
    TotalRealDemand::Float64
    TotalReactiveDemand::Float64

    # Damage data
    DISABLED::Vector{dataNode}
    HARDENED_DISABLED::Vector{dataNode}

    function problemData()
        p = new()

        p.numPhases = 3

        p.NODES = nodeData[]
        p.EDGES = edgeData[]
        p.GENERATORS = generatorData[]
        p.LOADS = loadData[]
        p.LINECODES = lineCodeData[]
        p.CYCLES = Vector{Int64}[]

        p.hashTableNodes = Dict{AbstractString, Int64}()
        p.hashTableEdges = Dict{AbstractString, Int64}()
        p.hashTableUniqueEdges = Dict{AbstractString, Int64}()
        p.hashTableGenerators = Dict{AbstractString, Int64}()
        p.hashTableLoads = Dict{AbstractString, Int64}()
        p.hashTableLineCodes = Dict{AbstractString, Int64}()

        p.LINE_CONSTRUCTION_COST = dataNode[]
        p.MICROGRID_COST = dataNode[]
        p.MICROGRID_FIXED_COST = dataNode[]
        p.MAX_MICROGRID = dataNode[]
        p.IS_CRITICAL_LOAD = dataNode[]
        p.HARDEN_COST = dataNode[]
        p.LINE_SWITCH_COST = dataNode[]

        p.CriticalRealDemand = 0.0
        p.CriticalReactiveDemand = 0.0
        p.TotalRealDemand = 0.0
        p.TotalReactiveDemand = 0.0

        p.DISABLED = dataNode[]
        p.HARDENED_DISABLED = dataNode[]

        loadProblemData(p)
        return p
    end
end

type ORDGDP     
    # MODEL
    mip_model::Model
    scen_idx::Int64

    # MASTER VARIABLES
    # TODO Get rid of the Any below!
    masterVariable::Vector{Any}

    # VARIABLES
    lineCycleVariable::Vector{Variable}
    switchCycleVariable::Vector{Variable}
    lineUseVariable::Vector{Variable}
    lineExistsVariable::Vector{Variable}
    lineDirectionVariableForward::Vector{Variable}
    lineDirectionVariableBackward::Vector{Variable}
    switchUseVariable::Vector{Variable}
    lineHardenVariable::Vector{Variable}
    flowRealVariable::Array{Variable,2}
    flowReactiveVariable::Array{Variable,2}
    voltageOffsetVariable::Array{Variable,2}
    voltageVariable::Array{Variable,2}
    loadRealVariable::Array{Variable,2}
    loadReactiveVariable::Array{Variable,2}
    generatorRealVariable::Array{Variable,2}
    generatorReactiveVariable::Array{Variable,2}
    loadServeVariable::Vector{Variable}
    facilityVariable::Vector{Variable}
   
    function ORDGDP(ordgdp_solver, problem_data::problemData, scen_idx::Int64)
        ordgdp = new()
        ordgdp.mip_model = Model(solver=ordgdp_solver)

        ordgdp.masterVariable = Any[]

        ordgdp.lineCycleVariable = Variable[]
        ordgdp.switchCycleVariable = Variable[]
        ordgdp.lineUseVariable = Variable[]
        ordgdp.lineExistsVariable = Variable[]
        ordgdp.lineDirectionVariableForward = Variable[]
        ordgdp.lineDirectionVariableBackward = Variable[]
        ordgdp.switchUseVariable = Variable[]
        ordgdp.lineHardenVariable = Variable[]
        ordgdp.flowRealVariable = Array{Variable,2}()
        ordgdp.flowReactiveVariable = Array{Variable,2}()
        ordgdp.voltageOffsetVariable = Array{Variable,2}()
        ordgdp.voltageVariable = Array{Variable,2}()
        ordgdp.loadRealVariable = Array{Variable,2}()
        ordgdp.loadReactiveVariable = Array{Variable,2}()
        ordgdp.generatorRealVariable = Array{Variable,2}()
        ordgdp.generatorReactiveVariable = Array{Variable,2}()
        ordgdp.loadServeVariable = Variable[]
        ordgdp.facilityVariable = Variable[]

        loadModel(ordgdp, problem_data, scen_idx)
        return ordgdp
    end
 
end


