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
    function nodeData(id::AbstractString, x::Float64, y::Float64, minVoltage::Float64, maxVoltage::Float64, refVoltage::Float64, hasPhase::Vector{Bool}, demand::Vector{Float64}, hasGenerator::Bool)
        n = new()
        n.id = id
        n.x = x
        n.y = y
        n.EdgeInList = Int64[] # derived from the edge data
        n.EdgeOutList = Int64[] # derived from the edge data
        n.GeneratorList = Int64[] # derived from the generator data
        n.LoadList = Int64[] # derived from the load data
        n.minVoltage = minVoltage
        n.maxVoltage = maxVoltage
        n.refVoltage = refVoltage
        n.hasPhase = hasPhase
        n.demand = demand # derived from the load data
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
    isNew::Bool
    
    function generatorData(id, node_id, hasPhase, maxRealPhase, maxReactivePhase, isNew)
        g = new()
        g.id = id
        g.node_id = node_id
        g.hasPhase = hasPhase
        g.maxRealPhase = maxRealPhase
        g.maxReactivePhase = maxReactivePhase
        g.isNew = isNew
        return g
    end
end

type loadData
    id::AbstractString
    node_id::AbstractString
    hasPhase::Vector{Bool}
    maxRealPhase::Vector{Float64}
    maxReactivePhase::Vector{Float64}
    function loadData(id::AbstractString, node_id::AbstractString, hasPhase::Vector{Bool}, maxRealPhase::Vector{Float64}, maxReactivePhase::Vector{Float64})
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
    lineCode::AbstractString
    numPoles::Int64
    isNew::Bool # whether or not this is a new line
    hasSwitch::Bool # whether or this line already has a switch 
    canHarden::Bool # whether or not we can harden this line
    canAddSwitch::Bool # whether or not we can add a switch here

    # Constructor    
    function edgeData(id, node1id, node2id, hasPhase, capacity, length, numPhases, isTransformer, lineCode, numPoles, isNew, hasSwitch, canHarden, canAddSwitch)
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
        e.isNew = isNew
        e.hasSwitch = hasSwitch
        e.canHarden = canHarden
        e.canAddSwitch = canAddSwitch
        return e
    end
end

type lineCodeData
    lineCode::AbstractString
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
    numScen::Int64

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
    DISABLED::Vector{Vector{dataNode}}
    HARDENED_DISABLED::Vector{Vector{dataNode}}

    # Constructor  TODO assumes the string is a filename. May want to make this abstract since it could just be a piped string
    function problemData(filename::AbstractString)
        p = new()

        p.numScen = 100

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

        p.DISABLED = Array(Vector{dataNode}, p.numScen)
        p.HARDENED_DISABLED = Array(Vector{dataNode}, p.numScen)
        for i in 1:p.numScen
            p.DISABLED[i] = dataNode[]
            p.HARDENED_DISABLED[i] = dataNode[]
        end

        loadProblemDataFile(p, filename)
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

# This function takes as input the problem data (which we will fill) and
# a filename of where the data comes from
function loadProblemDataFile(p::problemData, filename::AbstractString)
  data = JSON.parsefile(filename, dicttype=Dict, use_mmap=true)
  loadProblemDataJSONDict(p,data)
end

# This function takes as input the problem data (which we will fill) and
# a JSON string with the information we need
function loadProblemDataString(p::problemData, filename::AbstractString)
  data = JSON.parse(filename, dicttype = Dict{AbstractString,Any})
  loadProblemDataJSONDict(p,data)
end

# This function takes as input the problem data (which the function will fill) and
# a data dictionary with all the information we need
# TODO Have checks to see if data exists, and fill in defaults if it is not included... see the schema for the defaults
function loadProblemDataJSONDict(p::problemData, data::Dict)

  # load the the bus level data
  buses = data["buses"]
  for bus in buses
    id = bus["id"]
    x = bus["x"]  # TODO Emre is this field necessary?
    y = bus["y"] # TODO Emre is this field necessary?
    min_voltage = bus["min_voltage"]
    max_voltage = bus["max_voltage"]
    ref_voltage = bus["ref_voltage"]
    has_phase = [bus["has_phase"][1],bus["has_phase"][2],bus["has_phase"][3]] 
    demand = [0.0,0.0,0.0] # TODO Emre, not sure about this field?  Demand is real and reactive.      
    has_generator =  bus["has_generator"]
    node = nodeData(id, x, y, min_voltage, max_voltage, ref_voltage, has_phase, demand, has_generator)
    push!(p.NODES,node) 
    p.hashTableNodes[id] = length(p.NODES)
  end
    
  # grab the load data
  loads = data["loads"]
  for load in loads
    id = load["id"]
    node_id = load["node_id"]  
    has_phase = [load["has_phase"][1],load["has_phase"][2],load["has_phase"][3]]
    max_real_phase =  [load["max_real_phase"][1],load["max_real_phase"][2],load["max_real_phase"][3]]    
    max_reactive_phase =  [load["max_reactive_phase"][1],load["max_reactive_phase"][2],load["max_reactive_phase"][3]]   
    is_critical = false;
    if haskey(load, "is_critical")
      is_critical = load["is_critical"]
    end  
       
    l = loadData(id, node_id, has_phase, max_real_phase, max_reactive_phase)
    push!(p.LOADS,l) 
    p.hashTableLoads[id] = length(p.LOADS) 
    node_idx = p.hashTableNodes[node_id]
    push!(p.NODES[node_idx].LoadList, length(p.LOADS))    
    p.NODES[node_idx].demand = p.NODES[node_idx].demand + max_real_phase # TODO, is this what should be stored there?  

    push!(p.IS_CRITICAL_LOAD, dataNode(id, is_critical)) # TODO why not make this a member variable of edgeData?
  end

  # grab the generator data
  gens = data["generators"]
  for gen in gens
    id = gen["id"]
    node_id = gen["node_id"]
    has_phase = [gen["has_phase"][1],gen["has_phase"][2],gen["has_phase"][3]] 
    max_real_phase = [gen["max_real_phase"][1],gen["max_real_phase"][2],gen["max_real_phase"][3]]    
    max_reactive_phase = [gen["max_reactive_phase"][1],gen["max_reactive_phase"][2],gen["max_reactive_phase"][3]]           
    
    microgrid_cost = Inf
    if (haskey(gen, "microgrid_cost"))
      microgrid_cost = gen["microgrid_cost"]
    end
    
    microgrid_fixed_cost = Inf
    if (haskey(gen, "microgrid_fixed_cost"))
      microgrid_fixed_cost = gen["microgrid_fixed_cost"]
    end
    
    max_microgrid = 0
    if (haskey(gen, "max_micro_grid"))
      max_microgrid = gen["max_micro_grid"]
    end
    
    is_new = false    
    if (haskey(gen, "is_new"))
      is_new = gen["is_new"]
    end
      
      
    g = generatorData(id, node_id, has_phase, max_real_phase, max_reactive_phase, is_new)
    push!(p.GENERATORS,g) 
    p.hashTableGenerators[id] = length(p.GENERATORS) 
    node_idx = p.hashTableNodes[node_id]
    push!(p.NODES[node_idx].GeneratorList, length(p.GENERATORS))    
      
    push!(p.MICROGRID_COST, dataNode(id, microgrid_cost)) # TODO why not make this a member variable of edgeData?
    push!(p.MICROGRID_FIXED_COST, dataNode(id, microgrid_fixed_cost)) # TODO why not make this a member variable of edgeData?
    push!(p.MAX_MICROGRID, dataNode(id, max_microgrid)) # TODO why not make this a member variable of edgeData?
           
  end 
  
  # grab the edge data
  edges = data["lines"]
  for edge in edges
    id = edge["id"]
    node1id = edge["node1_id"]
    node2id = edge["node2_id"]
    has_phase = [edge["has_phase"][1],edge["has_phase"][2],edge["has_phase"][3]] 
    capacity = edge["capacity"]
    len = edge["length"]
    num_phases = edge["num_phases"]
    is_transformer = edge["is_transformer"]
    line_code = edge["line_code"]
    num_poles = edge["num_poles"] # TODO is this needed?
    construction_cost = edge["construction_cost"]
    
    harden_cost = Inf
    if (haskey(edge, "harden_cost"))
      harden_cost = edge["harden_cost"]         
    end
 
    switch_cost = Inf
    if (haskey(edge, "switch_cost"))
      harden_cost = edge["switch_cost"]         
    end
    
    is_new = false
    if (haskey(edge, "is_new"))
      is_new = edge["is_new"]         
    end
      
    has_switch = false
    if (haskey(edge, "has_switch"))
      has_switch = edge["has_switch"]         
    end
    
    can_harden = false
    if (haskey(edge, "can_harden"))
      can_harden = edge["can_harden"]         
    end
    
    can_add_switch = false
    if (haskey(edge, "can_add_switch"))
      can_add_switch = edge["can_add_switch"]         
    end
              
    e = edgeData(id, node1id, node2id, has_phase, capacity, len, num_phases, is_transformer, line_code, num_poles, is_new, has_switch, can_harden, can_add_switch)
  
    push!(p.EDGES,e) 
    p.hashTableEdges[id] = length(p.EDGES) 
    node1_idx = p.hashTableNodes[node1id]
    node2_idx = p.hashTableNodes[node2id]
    push!(p.NODES[node1_idx].EdgeInList, length(p.EDGES))    
    push!(p.NODES[node2_idx].EdgeOutList, length(p.EDGES))    
      
    push!(p.LINE_CONSTRUCTION_COST, dataNode(id, construction_cost)) # TODO why not make this a member variable of edgeData?
    push!(p.HARDEN_COST, dataNode(id, harden_cost)) # TODO why not make this a member variable of edgeData?
    push!(p.LINE_SWITCH_COST, dataNode(id, switch_cost)) # TODO why not make this a member variable of edgeData?
  end
  
  # grab the line code data
  line_codes = data["line_codes"]
  for line_code in line_codes
    id = line_code["line_code"]
    num_phases = line_code["num_phases"] # do we need this?  
    rmatrix = [ [line_code["rmatrix"][1][1], line_code["rmatrix"][1][2], line_code["rmatrix"][1][3]];
                [line_code["rmatrix"][2][1], line_code["rmatrix"][2][2], line_code["rmatrix"][2][3]];
                [line_code["rmatrix"][3][1], line_code["rmatrix"][3][2], line_code["rmatrix"][3][3]]
              ]

    xmatrix = [ [line_code["xmatrix"][1][1], line_code["xmatrix"][1][2], line_code["xmatrix"][1][3]];
                [line_code["xmatrix"][2][1], line_code["xmatrix"][2][2], line_code["xmatrix"][2][3]];
                [line_code["xmatrix"][3][1], line_code["xmatrix"][3][2], line_code["xmatrix"][3][3]]
              ]
    lc = lineCodeData(id, num_phases, rmatrix, xmatrix)
    push!(p.LINECODES,lc) 
    p.hashTableLineCodes[id] = length(p.LINECODES)     
  end  
    
  
       
    
  # Still need to hook up these things
  #   CYCLES  TODO how/when should this be populated?
  #   hashTableUniqueEdges::Dict{AbstractString,Int64}  TODO, how/when should this be populated
    
    # Demand data
  #  CriticalRealDemand::Float64
  #  CriticalReactiveDemand::Float64
  #  TotalRealDemand::Float64
  #  TotalReactiveDemand::Float64

  # Need to map (or add) the fields below to the fields above
  
  critical_load_met = data["critical_load_met"] # percentage of critical load met
  total_load_met = data["total_load_met"] # percenate of total load met
  chance_constraint = data["chance_constraint"] # chance constraint, for when we choose to add that back in
  phase_variation = data["phase_variation"] # for the phase variation at transformer constraint
        
  
    # Damage data TODO this does not seem the be the right data structures?  Vector of vectors?
  probDamage = 0.2
  for k in 1:p.numScen
    for i in 1:length(p.EDGES)
      for j in p.EDGES[i].numPoles
        dice = rand()
        if dice < probDamage
          push!(DISABLED[k], dataNode(p.EDGES[i].id, true))
        else
          push!(DISABLED[k], dataNode(p.EDGES[i].id, false))
       end 
       push!(HARDENED_DISABLED[k], dataNode(p.EDGES[i].id, false))
     end
  end

  # DETECT CYCLES
  G = Graph()
  for i in 1:length(p.NODES)
    addVertex(G, i)
  end
  for i in 1:length(p.EDGES)
    idx1 = p.hashTableVertex[p.EDGES[i].node1id]
    idx2 = p.hashTableVertex[p.EDGES[i].node2id]
    # Undirected graph
    addEdge(G, idx1, idx2)
    addEdge(G, idx2, idx1)
  end
  p.CYCLES = OrderedSet{Int64}[]

  detectCycles(G, p.CYCLES)

end
