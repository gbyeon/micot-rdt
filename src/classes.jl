# // -- ORDGDP CLASSES ---------------------------------

include("detect_cycles.jl")

# A simple struct for associating data with an identifier
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

# A struct containing node level information
type nodeData
    id::AbstractString
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
    function nodeData(id::AbstractString, minVoltage::Float64, maxVoltage::Float64, refVoltage::Float64, hasPhase::Vector{Bool}, demand::Vector{Float64}, hasGenerator::Bool)
        n = new()
        n.id = id
        n.EdgeInList = Int64[] # derived from the edge data
        n.EdgeOutList = Int64[] # derived from the edge data
        n.GeneratorList = Int64[] # derived from the generator data
        n.LoadList = Int64[] # derived from the load data
        n.minVoltage = minVoltage
        n.maxVoltage = maxVoltage
        n.refVoltage = refVoltage
        n.hasPhase = hasPhase
        n.demand = demand # derived from the load data, Emre what should go here
        n.hasGenerator = hasGenerator 
        return n
    end
end

# A struct containing information about generator data
type generatorData
    id::AbstractString
    node_id::AbstractString
    hasPhase::Vector{Bool}
    maxRealPhase::Vector{Float64}
    maxReactivePhase::Vector{Float64}
    isNew::Bool
    
  function generatorData(id::AbstractString, node_id::AbstractString, hasPhase::Vector{Bool}, maxRealPhase::Vector{Float64}, maxReactivePhase::Vector{Float64}, isNew::Bool)
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

#A struct containing information about load level data
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

# A struct containing information about edge level data
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
    isNew::Bool # whether or not this is a new line
    hasSwitch::Bool # whether or this line already has a switch 
    canHarden::Bool # whether or not we can harden this line
    canAddSwitch::Bool # whether or not we can add a switch here

    # Constructor    
    function edgeData(id::AbstractString, node1id::AbstractString, node2id::AbstractString, hasPhase::Vector{Bool}, capacity::Float64, length::Float64, numPhases::Int64, isTransformer::Bool, lineCode::AbstractString, isNew::Bool, hasSwitch::Bool, canHarden::Bool, canAddSwitch::Bool)
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
        e.isNew = isNew
        e.hasSwitch = hasSwitch
        e.canHarden = canHarden
        e.canAddSwitch = canAddSwitch
        return e
    end
end

# A struct containing information abouit line codes
type lineCodeData
    lineCode::AbstractString
    numPhases::Int64
    rmatrix::Array{Float64,3}
    xmatrix::Array{Float64,3}
    function lineCodeData(lineCode::AbstractString, numPhases::Int64, rmatrix::Array{Float64,3}, xmatrix::Array{Float64,3})
        l = new()
        l.lineCode = lineCode
        l.numPhases = numPhases
        l.rmatrix = rmatrix
        l.xmatrix = xmatrix
        return l
    end 
end

# a simple structure for scenario data
type scenarioData
  id::AbstractString
  disabled_edges::Dict{AbstractString,Bool}
  harden_disabled_edges::Dict{AbstractString,Bool}

  function scenarioData(id::AbstractString)
    s = new()
    s.id = id
    s.disabled_edges = Dict{AbstractString, Bool}()
    s.harden_disabled_edges = Dict{AbstractString, Bool}()    
    return s
  end
end
  
# A struct containing all the problem information
type problemData
    # DATA
    phaseVariation::Float64
    chanceConstraint::Float64

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

    # Scenario data
    SCENARIOS::Vector{scenarioData}
    
    # Constructor  RBENT TODO assumes the string is a filename. 
    # May want to make this abstract since it could just be a piped string
    function problemData(filename::AbstractString)
        p = new()

        p.numPhases = 3
        p.phaseVariation = .15;
        p.chanceConstraint = 1.0;

        p.NODES = nodeData[]
        p.EDGES = edgeData[]
        p.GENERATORS = generatorData[]
        p.LOADS = loadData[]
        p.LINECODES = lineCodeData[]
        p.SCENARIOS = scenarioData[]  
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

    # Constructor 
    # ordgdp_solver = optimization solver used to solve the problem
    # problem_data = the data of the power system
    # scen_idx = the scenario for which we are constructing the problem for
    # load_master = whether or not we are loading the master problem (base scenario) or not   
    function ORDGDP(ordgdp_solver, problem_data::problemData, scen_idx::Int64, load_master::Bool)
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

        loadModel(ordgdp, problem_data, scen_idx, load_master)
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
function loadProblemDataJSONDict(p::problemData, data::Dict)

  # load the the bus level data
  buses = data["buses"]
  for bus in buses
    id = bus["id"]
    min_voltage = haskey(bus,"min_voltage") ? bus["min_voltage"] : 0.8
    max_voltage = haskey(bus,"max_voltage") ? bus["max_voltage"] : 1.2
    ref_voltage = haskey(bus,"ref_voltage") ? bus["ref_voltage"] : 1.0
    has_phase = haskey(bus,"has_phase") ? [bus["has_phase"][1],bus["has_phase"][2],bus["has_phase"][3]] : [true, true, true]
    demand = [0.0,0.0,0.0] # TODO Emre, not sure about this field?  Demand is real and reactive.      
    has_generator =  haskey(bus,"has_generator") ? bus["has_generator"] : false
    node = nodeData(id, min_voltage, max_voltage, ref_voltage, has_phase, demand, has_generator)
    push!(p.NODES,node) 
    p.hashTableNodes[id] = length(p.NODES)
  end
  
  totalReal = 0
  totalReactive = 0
  totalCriticalReal = 0
  totalCriticalReactive = 0
    
  # grab the load data
  loads = data["loads"]
  for load in loads
    id = load["id"]
    node_id = load["node_id"]  
    has_phase = [load["has_phase"][1],load["has_phase"][2],load["has_phase"][3]]
    max_real_phase =  [load["max_real_phase"][1],load["max_real_phase"][2],load["max_real_phase"][3]]    
    max_reactive_phase =  [load["max_reactive_phase"][1],load["max_reactive_phase"][2],load["max_reactive_phase"][3]]   
    is_critical = haskey(load, "is_critical") ? load["is_critical"] : false;
       
    l = loadData(id, node_id, has_phase, max_real_phase, max_reactive_phase)
    push!(p.LOADS,l) 
    p.hashTableLoads[id] = length(p.LOADS) 
    node_idx = p.hashTableNodes[node_id]
    push!(p.NODES[node_idx].LoadList, length(p.LOADS))    
    p.NODES[node_idx].demand = p.NODES[node_idx].demand + max_real_phase # TODO, is this what should be stored there?  

    push!(p.IS_CRITICAL_LOAD, dataNode(id, is_critical)) # TODO why not make this a member variable of edgeData?

    if is_critical
      totalCriticalReal = totalCriticalReal + max_real_phase[1] + max_real_phase[2] + max_real_phase[3] 
      totalCriticalRective = totalCriticalReactive + max_reactive_phase[1] + max_reactive_phase[2] + max_reactive_phase[3]      
   else 
      totalReal = totalReal + max_real_phase[1] + max_real_phase[2] + max_real_phase[3] 
      totalRective = totalReactive + max_reactive_phase[1] + max_reactive_phase[2] + max_reactive_phase[3]          
    end
    
  end

  # grab the generator data
  gens = data["generators"]
  for gen in gens
    id = gen["id"]
    node_id = gen["node_id"]
    has_phase = haskey(gen,"has_phase") ? [gen["has_phase"][1],gen["has_phase"][2],gen["has_phase"][3]] : [true, true, true] 
    max_real_phase = haskey(gen,"max_real_phase") ? [gen["max_real_phase"][1],gen["max_real_phase"][2],gen["max_real_phase"][3]] : [Inf, Inf, Inf]   
    max_reactive_phase = haskey(gen,"max_reactive_phase") ? [gen["max_reactive_phase"][1],gen["max_reactive_phase"][2],gen["max_reactive_phase"][3]] : [Inf, Inf, Inf]              
    microgrid_cost = haskey(gen, "microgrid_cost") ? gen["microgrid_cost"] : Inf    
    microgrid_fixed_cost = haskey(gen, "microgrid_fixed_cost") ? gen["microgrid_fixed_cost"] : Inf    
    max_microgrid = haskey(gen, "max_micro_grid") ? gen["max_micro_grid"] : 0    
    is_new = haskey(gen, "is_new") ? gen["is_new"] : false    
            
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
    capacity = haskey(edge, "capacity") ? edge["capacity"] : Inf
    len = haskey(edge, "length") ? edge["length"] : 1.0
    num_phases = haskey(edge, "num_phases") ? edge["num_phases"] : 3
    is_transformer = haskey(edge, "is_transformer") ? edge["is_transformer"] : false
    line_code = edge["line_code"]
    is_new = haskey(edge, "is_new") ? edge["is_new"] : false        
    can_harden = haskey(edge, "can_harden") ? edge["can_harden"] : false    
    has_switch = haskey(edge, "has_switch") ? edge["has_switch"] : false    
    can_add_switch = haskey(edge, "can_add_switch") ? edge["can_add_switch"] : false
       
    construction_cost = haskey(edge,"construction_cost") ? edge["construction_cost"] : isNew ? Inf : 0.0    
    harden_cost = haskey(edge, "harden_cost") ? edge["harden_cost"] : can_harden ? 0 : Inf 
    switch_cost = haskey(edge, "switch_cost") ? edge["switch_cost"] : has_switch ? 0 : can_add_switch ? 0 : Inf    
              
    e = edgeData(id, node1id, node2id, has_phase, capacity, len, num_phases, is_transformer, line_code, is_new, has_switch, can_harden, can_add_switch)
  
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
            
#    rmatrix = [ line_code["rmatrix"][1][1] line_code["rmatrix"][1][2] line_code["rmatrix"][1][3];
 #               line_code["rmatrix"][2][1] line_code["rmatrix"][2][2] line_code["rmatrix"][2][3];
  #              line_code["rmatrix"][3][1] line_code["rmatrix"][3][2] line_code["rmatrix"][3][3]
   #           ]

#    xmatrix = [ line_code["xmatrix"][1][1] line_code["xmatrix"][1][2] line_code["xmatrix"][1][3];
 #               line_code["xmatrix"][2][1] line_code["xmatrix"][2][2] line_code["xmatrix"][2][3];
  #              line_code["xmatrix"][3][1] line_code["xmatrix"][3][2] line_code["xmatrix"][3][3]
   #           ]
      
      
    rmatrix = zeros(Float64, 3,3,3)  
    xmatrix = zeros(Float64, 3,3,3)
    for i=1:3
      for j=1:3
        rmatrix[1,i,j] = line_code["rmatrix"][i][j]          
        xmatrix[1,i,j] = line_code["xmatrix"][i][j]                   
      end
    end   
    
    # Rotate and scale the r and x matrices         
    # Index 1 = 120 and Index 2 = 240
    rotation = zeros(Float64, 2, 2, 2)        
    rotation[1,1,1] = -0.5;
    rotation[1,1,2] = -0.866;
    rotation[1,2,1] = 0.866;
    rotation[1,2,2] = -0.5;
    rotation[2,1,1] = -0.5;
    rotation[2,1,2] = 0.866;
    rotation[2,2,1] = -0.866;
    rotation[2,2,2] = -0.5;		

    for k=1:2
      for i=1:3
	for j=1:3
	  r = rmatrix[1,i,j]
	  x = xmatrix[1,i,j]
	  newr = (rotation[k,1,1] * r) + (rotation[k,1,2] * x) 
	  newx = (rotation[k,2,1] * r) + (rotation[k,2,2] * x) 
	  rmatrix[1+k,i,j] = newr
	  xmatrix[1+k,i,j] = newx
        end
     end
   end   
      
    lc = lineCodeData(id, num_phases, rmatrix, xmatrix)
    push!(p.LINECODES,lc) 
    p.hashTableLineCodes[id] = length(p.LINECODES)     
  end  
        
  # percentage of critical load met
  critical_load_met = haskey(data,"critical_load_met") ? data["critical_load_met"]  : 0.98
  
  # percentage of total load met
  total_load_met = haskey(data,"total_load_met") ? data["total_load_met"] : 0.5
    
  p.CriticalRealDemand = totalCriticalReal * critical_load_met    
  p.CriticalReactiveDemand = totalCriticalReactive * critical_load_met    
  p.TotalRealDemand = totalReal * total_load_met    
  p.TotalReactiveDemand = totalReactive * total_load_met      
  p.chanceConstraint = haskey(data,"chance_constraint") ? data["chance_constraint"] : 1.0
  p.phaseVariation = haskey(data,"phase_variation") ? data["phase_variation"] : 0.15
        
  # Grab the scenarios
  scenarios = data["scenarios"]
  for scenario in scenarios
    id = scenario["id"]
    sd = scenarioData(id)
    disabled_lines = scenario["disabled_lines"]
    harden_disabled_lines = scenario["hardened_disabled_lines"]  
    
    for disabled in disabled_lines
      sd.disabled_edges[disabled] = true
    end  
   
    for harden_disabled in harden_disabled_lines
      sd.harden_disabled_edges[harden_disabled] = true      
    end  
          
    push!(p.SCENARIOS,sd)
  end
    
  # DETECT CYCLES
  G = Graph()
  for i in 1:length(p.NODES)
    addVertex(G, i)
  end
  for i in 1:length(p.EDGES)
    idx1 = p.hashTableNodes[p.EDGES[i].node1id]
    idx2 = p.hashTableNodes[p.EDGES[i].node2id]
    # Undirected graph
    addEdge(G, idx1, idx2)
    addEdge(G, idx2, idx1)
  end
  p.CYCLES = OrderedSet{Int64}[]
  detectCycles(G, p.CYCLES)
  
  # index all edges by a the nodes that are connected to, uniquely
  unique_edge_idx = 1
  edge_assigned = Dict{Int64, Dict{Int64, Int64}}()
  for i in 1:length(p.NODES)
    edge_assigned[i] = Dict{Int64, Int64}()
  end
  for i in 1:length(p.EDGES)
    idx1 = p.hashTableNodes[p.EDGES[i].node1id]
    idx2 = p.hashTableNodes[p.EDGES[i].node2id]
    
    ei = -1
    if haskey(edge_assigned[idx1], idx2)
      ei = edge_assigned[idx1][idx2]
    end
    
    if ei == -1
      ei = unique_edge_idx
      unique_edge_idx = unique_edge_idx + 1
    end
    
    p.hashTableUniqueEdges[p.EDGES[i].id] = ei
    
  end  
    
  
end
