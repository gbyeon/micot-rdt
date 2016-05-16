// AUTHOR: EMRE YAMANGIL, 2014
// Contact: emreyamangil@gmail.com, emreyamangil@lanl.gov
// Description:
// 	Following code performs various algorithms to solve the Optimal Resilient Distribution Grid Design problem where every scenario corresponds to a subset of edges that are inoperable:
// 	1. Benders Scenario Based Decomposition,
// 	2. Variable Neighborhood Search,
// 	3. Benders Variable Neighborhood Decomposition Search.
// Note about CPLEX:
// 	The code heavily relies on CPLEX C++ API. For reference please consult ILOG CPLEX C++ API Reference Manual.

#include <ilcplex/ilocplex.h>
#include <map>
#include "read_text_data_damage.hpp"
#include "read_xml_data_damage.hpp"
#include "detect_cycles.hpp"
#include <numeric>
#include <omp.h>
#include <cmath>
#include <string>
//#include <algorithm>
//#include <Python.h>

ILOSTLBEGIN

#define RC_EPS 1.0e-6

struct Solver {
	IloEnv env;
	IloCplex cplex;
	Solver() : env(), cplex(env) {}
	void end() { env.end(); }
};

bool ismember(const vector<int> &vec, int x, int y) {
	bool ret_1 = false;
	bool ret_2 = false;
	for (int i = 0; i < vec.size(); ++i) {
		if (vec[i] == x) {
			ret_1 = true;
		} else if (vec[i] == y) {
			ret_2 = true;
		}
	}
	return ret_1 && ret_2;
}

template<class T> T vectorsum(const vector<T> &vec) {
	T sum = 0.0;
	for (int i = 0; i < vec.size(); ++i) {
		sum += vec[i];
	}
	return sum;
}

struct copyEdge {
	copyEdge() {}
	vector<string> id;
};

void addCycleConstraints(IloEnv env, vector<bool> include, IloRangeArray &cycleConstraint, IloNumVarArray &lineUseVariable, IloNumVarArray &switchUseVariable, const vector<int> &cycle, map<string, int> &hashTableEdge, map<string, int> &hashTableVertex, const vector<edgeData> &EDGES, int &numCon) {
	++numCon;
	//cout << numCon << endl;
	cycleConstraint.add(IloRange(env, -IloInfinity, (double) cycle.size() - 1.0));
	for (int k = 0; k < EDGES.size(); ++k) {
		if (ismember(cycle, hashTableVertex[EDGES[k].node1id], hashTableVertex[EDGES[k].node2id]) && include[k]) {
			cycleConstraint[numCon].setLinearCoef(lineUseVariable[k], 1.0);
			cycleConstraint[numCon].setLinearCoef(switchUseVariable[k], -1.0);
		} 
	}
}

void enumerateCycleConstraints(IloEnv env, vector<bool> include, int curIndex, IloRangeArray &cycleConstraint, IloNumVarArray &lineUseVariable, IloNumVarArray &switchUseVariable, const vector<int> &cycle, const vector<copyEdge> &CopyEdge, map<string, int> &hashTableEdge, map<string, int> &hashTableVertex, const vector<edgeData> &EDGES, int &numCon) {
	int ind = curIndex;
	for (;ind < include.size(); ++ind) {
		if (CopyEdge[ind].id.size() > 1) {
			break;
		} else {
			include[ind] = true;
		}
	}
	if (ind == include.size()) {
		addCycleConstraints(env, include, cycleConstraint, lineUseVariable, switchUseVariable, cycle, hashTableEdge, hashTableVertex, EDGES, numCon);
	} else {//if (CopyEdge[ind].id.size() > 1) {
		bool first = true;
		for (int i = 1; i < CopyEdge[ind].id.size(); ++i) {
			if (include[hashTableEdge[CopyEdge[ind].id[i]]]) {
				first = false;
			}
		}
		if (first) {
			include[ind] = true;
			enumerateCycleConstraints(env, include, ind+1, cycleConstraint, lineUseVariable, switchUseVariable, cycle, CopyEdge, hashTableEdge, hashTableVertex, EDGES, numCon);
			include[ind] = false;			
			enumerateCycleConstraints(env, include, ind+1, cycleConstraint, lineUseVariable, switchUseVariable, cycle, CopyEdge, hashTableEdge, hashTableVertex, EDGES, numCon);
		} else {
			include[ind] = false;			
			enumerateCycleConstraints(env, include, ind+1, cycleConstraint, lineUseVariable, switchUseVariable, cycle, CopyEdge, hashTableEdge, hashTableVertex, EDGES, numCon);
		}
	} /*else {
		include[ind] = true;
		enumerateCycleConstraints(env, include, ind+1, cycleConstraint, lineUseVariable, switchUseVariable, cycle, CopyEdge, hashTableEdge, hashTableVertex, EDGES, numCon);
	}*/
}

int minEdgeIndex (copyEdge &CopyEdge, map<string, int> &hashTableEdge) {
	int minInd = 10000;
	for (int i = 0; i < CopyEdge.id.size(); ++i) {
		if (minInd >= hashTableEdge[CopyEdge.id[i]]) {
			minInd = hashTableEdge[CopyEdge.id[i]];
		}
	}
	return minInd;
}

void GenerateObjectiveFunction(
	map<string, int> &hashTableVertex,
	map<string, int> &hashTableEdge,
	map<string, int> &hashTableGenerator,
	map<string, int> &hashTableLoad,
	vector<nodeData> &NODES,
	vector<edgeData> &EDGES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	vector<dataNode<double> > &LINE_CONSTRUCTION_COST,
	vector<dataNode<double> > &MICROGRID_COST,
	vector<dataNode<double> > &MICROGRID_FIXED_COST,
	vector<dataNode<bool> > &IS_CRITICAL_LOAD,
	vector<dataNode<double> > &MAX_MICROGRID,
	vector<dataNode<double> > &HARDEN_COST,
	vector<dataNode<double> > &LINE_SWITCH_COST,
	IloObjective &SubObjective,
	vector<double> &ObjectiveVector,
	IloNumVarArray &lineUseVariable, // VARIABLES
	IloNumVarArray &lineExistsVariable, // VARIABLES
	IloNumVarArray &lineDirectionVariable,
	IloNumVarArray &lineCycleVariable,
	IloNumVarArray &lineHardenVariable,
	IloNumVarArray &switchUseVariable,
	IloNumVarArray &switchCycleVariable,
	IloNumVarArray &flowVariable,
	IloNumVarArray &flowReactiveVariable,
	IloNumVarArray &voltageVariable,
	IloNumVarArray &voltageOffsetVariable,
	IloNumVarArray &loadVariable,
	IloNumVarArray &loadReactiveVariable,
	IloNumVarArray &loadServeVariable,
	IloNumVarArray &generatorVariable,
	IloNumVarArray &generatorReactiveVariable,
	//IloNumVarArray &microgridVariable,
	IloNumVarArray &facilityVariable) {


	vector<double> MAX_ADDED(GENERATORS.size(), 0.0);
	for (int j = 0; j < MAX_MICROGRID.size(); ++j) {
		MAX_ADDED[hashTableGenerator[MAX_MICROGRID[j].id]] = MAX_MICROGRID[j].data;
	}


	vector<double> FACILITY_COST(GENERATORS.size(), 0.0);
	for (int j = 0; j < MICROGRID_COST.size(); ++j) {
		FACILITY_COST[hashTableGenerator[MICROGRID_COST[j].id]] += MAX_ADDED[hashTableGenerator[MICROGRID_COST[j].id]] * MICROGRID_COST[j].data;
	}
	for (int j = 0; j < MICROGRID_FIXED_COST.size(); ++j) {
		FACILITY_COST[hashTableGenerator[MICROGRID_FIXED_COST[j].id]] += MICROGRID_FIXED_COST[j].data;
	}
	for (int j = 0; j < LINE_CONSTRUCTION_COST.size(); ++j) {
		ObjectiveVector[hashTableEdge[LINE_CONSTRUCTION_COST[j].id]] = 
			LINE_CONSTRUCTION_COST[j].data;
		SubObjective.setLinearCoef(lineUseVariable[hashTableEdge[LINE_CONSTRUCTION_COST[j].id]], 
			LINE_CONSTRUCTION_COST[j].data);
			//(double) round(LINE_CONSTRUCTION_COST[j].data * EDGES[j].length + 5.0));
	}

	for (int j = 0; j < HARDEN_COST.size(); ++j) {
		ObjectiveVector[2 * EDGES.size() + hashTableEdge[HARDEN_COST[j].id]] = HARDEN_COST[j].data;
		SubObjective.setLinearCoef(lineHardenVariable[hashTableEdge[HARDEN_COST[j].id]], HARDEN_COST[j].data);
	}
	/*for (int j = 0; j < EDGES.size(); ++j) {
		ObjectiveVector[1 * EDGES.size() + j] = 25.0;
		SubObjective.setLinearCoef(switchUseVariable[j], 25.0);
	}*/
	for (int j = 0; j < LINE_SWITCH_COST.size(); ++j) {
		ObjectiveVector[1 * EDGES.size() + hashTableEdge[LINE_SWITCH_COST[j].id]] = LINE_SWITCH_COST[j].data;
		SubObjective.setLinearCoef(switchUseVariable[hashTableEdge[LINE_SWITCH_COST[j].id]], LINE_SWITCH_COST[j].data);
	}
	for (int j = 0; j < GENERATORS.size(); ++j) {
		ObjectiveVector[3 * EDGES.size() + j] = FACILITY_COST[j];
		SubObjective.setLinearCoef(facilityVariable[j], FACILITY_COST[j]);
	}



}



// FOLLOWING METHOD POPULATES A SINGLE SCENARIO
void PopulateScenarioData(
	Solver &SubSolver,
	map<string, int> &hashTableVertex,
	map<string, int> &hashTableEdge,
	map<string, int> &hashTableGenerator,
	map<string, int> &hashTableLoad,
	vector<nodeData> &NODES,
	vector<edgeData> &EDGES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	map<int, lineCodeData> &LINECODES,
	vector<dataNode<double> > &LINE_CONSTRUCTION_COST,
	vector<dataNode<double> > &MICROGRID_COST,
	vector<dataNode<double> > &MICROGRID_FIXED_COST,
	vector<dataNode<bool> > &IS_CRITICAL_LOAD,
	vector<dataNode<double> > &MAX_MICROGRID,
	vector<dataNode<double> > &HARDEN_COST,
	vector<dataNode<double> > &LINE_SWITCH_COST,
	double LoadMet,
	double CriticalLoadMet,
	double PhaseVariation,
	vector<graph_vertex> &G,
	vector<vector<int> > &CYCLES,
	IloObjective &SubObjective,
	vector<double> &ObjectiveVector,
	IloRangeArray &lineExistsConstraint, // CONSTRAINTS
	IloRangeArray &lineConstraint1, // CONSTRAINTS
	IloRangeArray &lineConstraint2,
	IloRangeArray &lineReactiveConstraint1,
	IloRangeArray &lineReactiveConstraint2,
	IloRangeArray &directionConstraint,
	IloRangeArray &switchConstraint1,
	IloRangeArray &switchConstraint2,
	IloRangeArray &switchReactiveConstraint1,
	IloRangeArray &switchReactiveConstraint2,
	IloRangeArray &variationConstraint1,
	IloRangeArray &variationConstraint2,
	IloRangeArray &variationReactiveConstraint1,
	IloRangeArray &variationReactiveConstraint2,
	IloRangeArray &damageConstraint,
	IloRangeArray &loadConstraint,
	IloRangeArray &loadReactiveConstraint,
	IloRangeArray &generationConstraint,
	IloRangeArray &generationReactiveConstraint,
	IloRangeArray &balanceConstraint,
	IloRangeArray &balanceReactiveConstraint,
	IloRangeArray &microgridConstraint,
	IloRangeArray &cycleConstraint,
	IloRangeArray &cycleSubConstraint,
	IloRangeArray &switchableConstraint,
	IloRangeArray &criticalServeConstraint,
	IloRangeArray &criticalServeReactiveConstraint,
	IloRangeArray &loadServeConstraint,
	IloRangeArray &loadServeReactiveConstraint,
	IloRangeArray &switchImplicationConstraint1,
	//IloRangeArray &switchImplicationConstraint2,
	IloRangeArray &hardenConstraint,
	IloRangeArray &voltageConstraint,
	IloRangeArray &voltageOffsetConstraint,
	IloRangeArray &voltageVariationConstraint1,
	IloRangeArray &voltageVariationConstraint2,
	IloNumVarArray &lineUseVariable, // VARIABLES
	IloNumVarArray &lineExistsVariable, // VARIABLES
	IloNumVarArray &lineDirectionVariable,
	IloNumVarArray &lineCycleVariable,
	IloNumVarArray &lineHardenVariable,
	IloNumVarArray &switchUseVariable,
	IloNumVarArray &switchCycleVariable,
	IloNumVarArray &flowVariable,
	IloNumVarArray &flowReactiveVariable,
	IloNumVarArray &voltageVariable,
	IloNumVarArray &voltageOffsetVariable,
	IloNumVarArray &loadVariable,
	IloNumVarArray &loadReactiveVariable,
	IloNumVarArray &loadServeVariable,
	IloNumVarArray &generatorVariable,
	IloNumVarArray &generatorReactiveVariable,
	//IloNumVarArray &microgridVariable,
	IloNumVarArray &facilityVariable,
	int &NumDamagedEdges,
	int NumPhases,
	int ScenarioIndex,
	int NumScenarios,
	double discreteMicrogrid,
	string damage_str,
	string root_str,
	bool HardenedDisabled,
	vector<dataNode<bool> > &DISABLED,
	vector<dataNode<bool> > &HARDENED_DISABLED,
	double newImpedanceMultiplier,
	bool LDFIndicator) {
	cout << "// -- STARTED NEW SCENARIO #" << ScenarioIndex+1 << " ----------------------------------" << endl;

	// GENERATES THE SET OF NODES WHERE WE NEED TO WRITE FLOW BALANCE FOR A PARTICULAR PHASE
	vector<vector<bool> > balancable(NODES.size(), vector<bool>(NumPhases, false));
	for (int j = 0; j < NODES.size(); ++j) {
		for (list<string>::iterator it = G[j].EdgeID.begin(); it != G[j].EdgeID.end(); ++it) {
			for (int k = 0; k < NumPhases; ++k) {
				if (EDGES[hashTableEdge[*it]].hasphase[k]) {		
					balancable[j][k] = true;
				}
			}
		}
	}

	// FIND COPY EDGES
	vector<copyEdge> CopyEdge(EDGES.size());
	for (int j = 0; j < EDGES.size(); ++j) {
		CopyEdge[j].id.push_back(EDGES[j].id);
	}
	for (int j = 0; j < EDGES.size(); ++j) {
		//if (EDGES[j].istransformer) cout << EDGES[j].id << endl;
		for (int k = j+1; k < EDGES.size(); ++k) {
			if (hashTableVertex[EDGES[j].node1id] == hashTableVertex[EDGES[k].node1id] && hashTableVertex[EDGES[j].node2id] == hashTableVertex[EDGES[k].node2id]) {
				CopyEdge[j].id.push_back(EDGES[k].id);
				CopyEdge[k].id.push_back(EDGES[j].id);
			} 
		}
	}
	// SECOND WAY OF ENUMERATING CYCLES
	map<string, int> hashTableCopyEdge;
	int edgeInd = -1;
	for (int i = 0; i < EDGES.size(); ++i) {
		int minInd = minEdgeIndex(CopyEdge[i], hashTableEdge);
		if (minInd == i) {
			++edgeInd;
			hashTableCopyEdge.insert(pair<string, int>(EDGES[i].id, edgeInd));
		} else {
			hashTableCopyEdge.insert(pair<string, int>(EDGES[i].id, hashTableCopyEdge[EDGES[minEdgeIndex(CopyEdge[i], hashTableEdge)].id]));
		}
	}
	int NumDistinctEdges = edgeInd + 1;

	IloEnv env = SubSolver.cplex.getEnv();
	
	std::stringstream ss;
	ss << ScenarioIndex + 1;
	std::string str = ss.str();
	string input_file = root_str + damage_str + "/config-34bussnow-" + str + ".xml";

	// GENERATES THE SET OF EDGES WHICH MIGHT BE HARDENED
	vector<bool> supportable(EDGES.size(), false);
	for (int j = 0; j < HARDEN_COST.size(); ++j) {
		supportable[hashTableEdge[HARDEN_COST[j].id]] = true;
	}
        
	vector<bool> edgeExists(EDGES.size(), true);
	for (int j = 0; j < LINE_CONSTRUCTION_COST.size(); ++j) {
		if (LINE_CONSTRUCTION_COST[j].data != 0.0) {
			edgeExists[hashTableEdge[LINE_CONSTRUCTION_COST[j].id]] = false;
		}
	}

	vector<bool> switchable(EDGES.size(), false);
	for (int j = 0; j < LINE_SWITCH_COST.size(); ++j) {
		//if (LINE_SWITCH_COST[j].data >= 0.0) {
			switchable[hashTableEdge[LINE_SWITCH_COST[j].id]] = true;
		//}
	}
	// -- VARIABLES -------------------------------------------------------------------
	// CYCLE VARIABLES
	for (int j = 0; j < NumDistinctEdges; ++j) {
		std::stringstream ss1;
		ss1 << j;
		str = "lineCycleVariable_" + ss1.str();
		lineCycleVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str()));
		str = "switchCycleVariable_" + ss1.str();
		switchCycleVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str()));
	}
	// EDGE BASED VARIABLES
	// For now I added all possible variables (if it doesn't has phase, it adds a 0 constant, CPLEX PreSolve filters constants)
	for (int j = 0; j < EDGES.size(); ++j) {
		std::stringstream ss1;
		ss1 << EDGES[j].id;
		str = "lineUseVariable_" + ss1.str();
		lineUseVariable.add(IloNumVar(env, edgeExists[j] ? 1 : 0, 1, ILOINT, str.c_str())); // LB: edgeExists[j] ? 1 : 0
		//lineUseVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str())); // LB: edgeExists[j] ? 1 : 0
		str = "lineExistsVariable_" + ss1.str();
		lineExistsVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str())); // LB: edgeExists[j] ? 1 : 0
		str = "lineDirectionVariable_" + ss1.str() + "_0";
		lineDirectionVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str()));
		str = "lineDirectionVariable_" + ss1.str() + "_1";
		lineDirectionVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str()));
		str = "switchUseVariable_" + ss1.str();
		switchUseVariable.add(IloNumVar(env, 0, switchable[j] ? 1 : 0, ILOINT, str.c_str()));
		//switchUseVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str())); // ALLOW SWITCHES EVERYWHERE
		str = "lineHardenVariable_" + ss1.str();
		lineHardenVariable.add(IloNumVar(env, 0, supportable[j] ? 1 : 0, ILOINT, str.c_str())); // UB: supportable[j] ? 1 : 0
		//lineHardenVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str())); // UB: supportable[j] ? 1 : 0
		for (int k = 0; k < NumPhases; ++k) {
			std::stringstream ss2;
			ss2 << k;
			str = "flowVariable_" + ss1.str() + "_" + ss2.str();
			//cout << str << endl;
			flowVariable.add(IloNumVar(env, EDGES[j].hasphase[k] ? -IloInfinity : 0.0, EDGES[j].hasphase[k] ? IloInfinity : 0.0, ILOFLOAT, str.c_str()));
			str = "flowReactiveVariable_" + ss1.str() + "_" + ss2.str();
			flowReactiveVariable.add(IloNumVar(env, EDGES[j].hasphase[k] ? -IloInfinity : 0.0, EDGES[j].hasphase[k] ? IloInfinity : 0.0, ILOFLOAT, str.c_str()));
			str = "voltageOffsetVariable_" + ss1.str() + "_" + ss2.str();
			voltageOffsetVariable.add(IloNumVar(env, -IloInfinity, IloInfinity, ILOFLOAT, str.c_str())); // NODES[j].hasphase[k] ? IloInfinity : 0.0
		}
	}
	// NODE BASED VARIABLES
	// For every node
	for (int j = 0; j < NODES.size(); ++j) {
		std::stringstream ss1;
		ss1 << NODES[j].id;
		for (int k = 0; k < NumPhases; ++k) {
			std::stringstream ss2;
			ss2 << k;
			str = "voltageVariable_" + ss1.str() + "_" + ss2.str();
			voltageVariable.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT, str.c_str())); // NODES[j].hasphase[k] ? IloInfinity : 0.0
			str = "loadVariable_" + ss1.str() + "_" + ss2.str();
			loadVariable.add(IloNumVar(env, 0.0, NODES[j].hasphase[k] ? IloInfinity : 0.0, ILOFLOAT, str.c_str()));
			str = "loadReactiveVariable_" + ss1.str() + "_" + ss2.str();
			loadReactiveVariable.add(IloNumVar(env, 0.0, NODES[j].hasphase[k] ? IloInfinity : 0.0, ILOFLOAT, str.c_str()));
			str = "generatorVariable_" + ss1.str() + "_" + ss2.str();
			generatorVariable.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT, str.c_str()));		
			str = "generatorReactiveVariable_" + ss1.str() + "_" + ss2.str();
			generatorReactiveVariable.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT, str.c_str()));		
			//cout << (j * NumPhases + k) << " " << (NODES[j].hasphase[k] ? IloInfinity : 0.0) << endl;
		}
	}	
	// LOAD BASED VARIABLES
	// For every load
	for (int j = 0; j < LOADS.size(); ++j) {
		std::stringstream ss1;
		ss1 << LOADS[j].id;
		str = "loadServeVariable_" + ss1.str();
		loadServeVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str()));
	}
	// GENERATOR VARIABLES
	// For every generator
	for (int j = 0; j < GENERATORS.size(); ++j) {
		std::stringstream ss1;
		ss1 << GENERATORS[j].id;
		str = "facilityVariable_" + ss1.str();
		facilityVariable.add(IloNumVar(env, 0, 1, ILOINT, str.c_str()));
		/*for (int k = 0; k < NumPhases; ++k) {
			std::stringstream ss2;
			ss2 << k;
			str = "microgridVariable_" + ss1.str() + "_" + ss2.str();
			microgridVariable.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT, str.c_str()));
		}*/
	}
	// -- CONSTRAINTS -----------------------------------------------------------------
	// VOLTAGE CONSTRAINTS
	// RESISTANCE AND REACTANCE VALUES ARE MULTIPLIED BY 5.28 IN READ TEXT DATA
	// 163264 213007 221416
	//voltageVariable[hashTableVertex["sourcebus"]+0].setUB(pow(69.0, 2.0));
	//voltageVariable[hashTableVertex["sourcebus"]+1].setUB(pow(69.0, 2.0));
	//voltageVariable[hashTableVertex["sourcebus"]+2].setUB(pow(69.0, 2.0));
	//voltageVariable[hashTableVertex["sourcebus"]+0].setLB(pow(69.0, 2.0));
	//voltageVariable[hashTableVertex["sourcebus"]+1].setLB(pow(69.0, 2.0));
	//voltageVariable[hashTableVertex["sourcebus"]+2].setLB(pow(69.0, 2.0));
	int numCon = -1;
	if (LDFIndicator) {
		for (int j = 0; j < EDGES.size(); ++j) {
			if (!EDGES[j].istransformer) {
				bool multiplierCondi = (EDGES[j].id[0] == 'n' || EDGES[j].id[0] == 'o');//(hashTableEdge["l2001"] == j);
				double impedanceMultiplier = (multiplierCondi ? newImpedanceMultiplier : 1.0); //2000 works
				//double impedanceMultiplier = 1.0;
				if (EDGES[j].NumPhases == 1) {
					++numCon;
					voltageConstraint.add(IloRange(env, 0.0, 0.0));
					int phaseind = -1;
					for (int k = 0; k < NumPhases; ++k) {
						if (EDGES[j].hasphase[k]) {
							phaseind = k;
							break;
						}
					}
					voltageConstraint[numCon].setLinearCoef(voltageVariable[
						hashTableVertex[EDGES[j].node1id]*NumPhases + phaseind], 1.0);
					voltageConstraint[numCon].setLinearCoef(voltageVariable[
						hashTableVertex[EDGES[j].node2id]*NumPhases + phaseind], -1.0);
					voltageConstraint[numCon].setLinearCoef(voltageOffsetVariable[
						j*NumPhases + phaseind], 1.0);


					map<int, lineCodeData>::iterator it = LINECODES.find(EDGES[j].linecode);

					voltageConstraint[numCon].setLinearCoef(flowVariable[
						j*NumPhases + phaseind],
						-2.0 * it->second.rmatrix[0][0][0] * impedanceMultiplier); // MAGIC LINES
					voltageConstraint[numCon].setLinearCoef(flowReactiveVariable[
						j*NumPhases + phaseind],
						-2.0 * it->second.xmatrix[0][0][0] * impedanceMultiplier); // MAGIC LINES


				} else if (EDGES[j].NumPhases == 3) {
					for (int k = 0; k < NumPhases; ++k) {
						if (EDGES[j].hasphase[k]) {
							++numCon;
							voltageConstraint.add(IloRange(env, 0.0, 0.0));

							voltageConstraint[numCon].setLinearCoef(voltageVariable[
								hashTableVertex[EDGES[j].node1id]*NumPhases + k], 1.0);
							voltageConstraint[numCon].setLinearCoef(voltageVariable[
								hashTableVertex[EDGES[j].node2id]*NumPhases + k], -1.0);
							voltageConstraint[numCon].setLinearCoef(voltageOffsetVariable[
								j*NumPhases + k], 1.0);

							map<int, lineCodeData>::iterator it = LINECODES.find(EDGES[j].linecode);
							for (int l = 0; l < NumPhases; ++l) {
								if (EDGES[j].hasphase[l]) {
									voltageConstraint[numCon].setLinearCoef(flowVariable[
										j*NumPhases + l],
										-2.0 * it->second.rmatrix[(3-k+l)%NumPhases][k][l] * impedanceMultiplier); // MAGIC LINES

									voltageConstraint[numCon].setLinearCoef(flowReactiveVariable[
										j*NumPhases + l],
										-2.0 * it->second.xmatrix[(3-k+l)%NumPhases][k][l] * impedanceMultiplier); // MAGIC LINES

								}
							}
						}
					}
				}
			}
		}
	}
	numCon = -1;
	double voltageMAX = 1e+6;
	for (int j = 0; j < EDGES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			// voltageOffsetVarible <= voltageMAX (1 - (lineUseVariable - switchUseVariable)) = voltageMAX ( 1 - (lineExistsVariable));
			// voltageOffsetVarible >= -voltageMAX (1 - (lineUseVariable - switchUseVariable));
			++numCon;
			voltageOffsetConstraint.add(IloRange(env, -IloInfinity, voltageMAX));
			voltageOffsetConstraint[numCon].setLinearCoef(voltageOffsetVariable[j * NumPhases + k], 1.0);
			voltageOffsetConstraint[numCon].setLinearCoef(lineExistsVariable[j], voltageMAX);
			//voltageOffsetConstraint[numCon].setLinearCoef(switchUseVariable[j], -voltageMAX);
			++numCon;
			voltageOffsetConstraint.add(IloRange(env, -voltageMAX, IloInfinity));
			voltageOffsetConstraint[numCon].setLinearCoef(voltageOffsetVariable[j * NumPhases + k], 1.0);
			voltageOffsetConstraint[numCon].setLinearCoef(lineExistsVariable[j], -voltageMAX);
			//voltageOffsetConstraint[numCon].setLinearCoef(switchUseVariable[j], voltageMAX);
		}
	}
	numCon = -1;
	for (int j = 0; j < EDGES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			++numCon;
			lineExistsConstraint.add(IloRange(env, 0.0, 0.0));
			lineExistsConstraint[numCon].setLinearCoef(lineExistsVariable[j], 1.0);
			lineExistsConstraint[numCon].setLinearCoef(lineUseVariable[j], -1.0);
			lineExistsConstraint[numCon].setLinearCoef(switchUseVariable[j], 1.0);
		}
	}
	numCon = -1;
	for (int j = 0; j < NODES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			//if (NODES[j].hasphase[k]) {// || hashTableVertex["2800"] == j || hashTableVertex["1800"] == j || hashTableVertex["800"] == j) {
				// voltageVar <= voltageRef * 1.05^2 
				// voltageVar >= voltageRef * 0.95^2 
				++numCon;
				voltageVariationConstraint1.add(IloRange(env, 
					-IloInfinity, pow(NODES[j].refVoltage * 1.05, 2.0) * 1e+3)); // 1.05!!!!
				// * 1e+3 if we are dealing with kW
				if (j != 0) voltageVariationConstraint2.add(IloRange(env, 
					pow(NODES[j].refVoltage * 0.95, 2.0) * 1e+3, IloInfinity)); // 0.95!!!!
				else voltageVariationConstraint2.add(IloRange(env, 
					pow(NODES[j].refVoltage * 1.05, 2.0) * 1e+3, IloInfinity)); // 0.95!!!!
				// * 1e+3 if we are dealing with kW

				voltageVariationConstraint1[numCon].setLinearCoef(
					voltageVariable[j*NumPhases + k], 1.0);

				voltageVariationConstraint2[numCon].setLinearCoef(
					voltageVariable[j*NumPhases + k], 1.0);

				/*voltageVariationConstraint1.add(IloRange(env, -IloInfinity, 0.0));
				voltageVariationConstraint2.add(IloRange(env, -IloInfinity, 0.0));
				
				voltageVariationConstraint1[numCon].setLinearCoef(
					voltageVariable[j*NumPhases + k], 1.0);
				voltageVariationConstraint1[numCon].setLinearCoef(
					voltageVariable[hashTableVertex["sourcebus"] * NumPhases + k], - pow(1.05, 2.0));

				voltageVariationConstraint2[numCon].setLinearCoef(
					voltageVariable[j*NumPhases + k], -1.0);
				voltageVariationConstraint2[numCon].setLinearCoef(
					voltageVariable[hashTableVertex["sourcebus"] * NumPhases + k], pow(0.95, 2.0));*/

			//}
		}
	}

	// LINE CONSTRAINTS
	// FROM EDGES[j].capacity * (1.0 + PhaseVariation) TO EDGES[j].capacity / EDGES[j].NumPhases
	numCon = -1;
	for (int j = 0; j < EDGES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			if (EDGES[j].hasphase[k]) {
				++numCon;
				lineConstraint1.add(IloRange(env, 0.0, IloInfinity));
				lineConstraint2.add(IloRange(env, 0.0, IloInfinity));
				lineConstraint1[numCon].setLinearCoef(flowVariable[j*NumPhases + k], 1.0);
				lineConstraint2[numCon].setLinearCoef(flowVariable[j*NumPhases + k], -1.0);
				lineConstraint1[numCon].setLinearCoef(lineDirectionVariable[2*j], EDGES[j].capacity / EDGES[j].NumPhases);
				lineConstraint2[numCon].setLinearCoef(lineDirectionVariable[2*j+1], EDGES[j].capacity / EDGES[j].NumPhases);		
			}
		}		
	}
	numCon = -1;
	for (int j = 0; j < EDGES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			if (EDGES[j].hasphase[k]) {
				++numCon;
				lineReactiveConstraint1.add(IloRange(env, 0.0, IloInfinity));
				lineReactiveConstraint2.add(IloRange(env, 0.0, IloInfinity));
				lineReactiveConstraint1[numCon].setLinearCoef(flowReactiveVariable[j*NumPhases + k], 1.0);
				lineReactiveConstraint2[numCon].setLinearCoef(flowReactiveVariable[j*NumPhases + k], -1.0);
				lineReactiveConstraint1[numCon].setLinearCoef(lineDirectionVariable[2*j], EDGES[j].capacity / EDGES[j].NumPhases);
				lineReactiveConstraint2[numCon].setLinearCoef(lineDirectionVariable[2*j+1], EDGES[j].capacity / EDGES[j].NumPhases);		
			}
		}		
	}
	for (int j = 0; j < EDGES.size(); ++j) {
		directionConstraint.add(IloRange(env, -IloInfinity, 0.0));
		directionConstraint[j].setLinearCoef(lineExistsVariable[j], -1.0); // CHANGED FROM lineUseVariable
		directionConstraint[j].setLinearCoef(lineDirectionVariable[2*j], 1.0);
		directionConstraint[j].setLinearCoef(lineDirectionVariable[2*j+1], 1.0);
	}
	cout << "Number of Distinct Edges: " << NumDistinctEdges << endl;
	// SWITCH CONSTRAINTS
	// FROM EDGES[j].capacity * (1.0 + PhaseVariation) TO EDGES[j].capacity / EDGES[j].NumPhases
	/*numCon = -1;
	for (int j = 0; j < EDGES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			if (EDGES[j].hasphase[k]) {
				++numCon;
				switchConstraint1.add(IloRange(env, -EDGES[j].capacity / EDGES[j].NumPhases, IloInfinity));
				switchConstraint2.add(IloRange(env, -EDGES[j].capacity / EDGES[j].NumPhases, IloInfinity));
				switchConstraint1[numCon].setLinearCoef(flowVariable[j*NumPhases + k], 1.0);
				switchConstraint2[numCon].setLinearCoef(flowVariable[j*NumPhases + k], -1.0);
				switchConstraint1[numCon].setLinearCoef(switchUseVariable[j], -EDGES[j].capacity / EDGES[j].NumPhases);
				switchConstraint2[numCon].setLinearCoef(switchUseVariable[j], -EDGES[j].capacity / EDGES[j].NumPhases);		
			}
		}		
	}
	numCon = -1;
	for (int j = 0; j < EDGES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			if (EDGES[j].hasphase[k]) {
				++numCon;
				switchReactiveConstraint1.add(IloRange(env, -EDGES[j].capacity / EDGES[j].NumPhases, IloInfinity));
				switchReactiveConstraint2.add(IloRange(env, -EDGES[j].capacity / EDGES[j].NumPhases, IloInfinity));
				switchReactiveConstraint1[numCon].setLinearCoef(flowReactiveVariable[j*NumPhases + k], 1.0);
				switchReactiveConstraint2[numCon].setLinearCoef(flowReactiveVariable[j*NumPhases + k], -1.0);
				switchReactiveConstraint1[numCon].setLinearCoef(switchUseVariable[j], -EDGES[j].capacity / EDGES[j].NumPhases);
				switchReactiveConstraint2[numCon].setLinearCoef(switchUseVariable[j], -EDGES[j].capacity / EDGES[j].NumPhases);		
			}
		}		
	}*/
	// VARIATION CONSTRAINTS
	//PhaseVariation = 0.5;
	numCon = -1;
	for (int j = 0; j < EDGES.size(); ++j) {
		if (EDGES[j].NumPhases > 1 && EDGES[j].istransformer) { // No need for single phase flow
			if (EDGES[j].NumPhases == 2) cout << "EDGE WITH 2 PHASE FLOW!" << endl;
			for (int k = 0; k < NumPhases; ++k) {
				if (EDGES[j].hasphase[k]) {
					++numCon;
					variationConstraint1.add(IloRange(env, 0.0, IloInfinity));
					variationConstraint2.add(IloRange(env, 0.0, IloInfinity));
					for (int l = 0; l < NumPhases; ++l) {
						if (k != l && EDGES[j].hasphase[l]) {
							variationConstraint1[numCon].setLinearCoef(flowVariable[j*NumPhases + l], (PhaseVariation - 1.0) / (double) EDGES[j].NumPhases);
							variationConstraint2[numCon].setLinearCoef(flowVariable[j*NumPhases + l], (PhaseVariation + 1.0) / (double) EDGES[j].NumPhases);
						}
					}	
					variationConstraint1[numCon].setLinearCoef(flowVariable[j*NumPhases + k], ((double) EDGES[j].NumPhases + PhaseVariation - 1.0) / (double) EDGES[j].NumPhases);
					variationConstraint2[numCon].setLinearCoef(flowVariable[j*NumPhases + k], (PhaseVariation + 1.0 - (double) EDGES[j].NumPhases) / (double) EDGES[j].NumPhases);		
				}
			}
		}
	}
	numCon = -1;
	for (int j = 0; j < EDGES.size(); ++j) {
		if (EDGES[j].NumPhases > 1 && EDGES[j].istransformer) { // No need for single phase flow
			if (EDGES[j].NumPhases == 2) cout << "EDGE WITH 2 PHASE FLOW!" << endl;
			for (int k = 0; k < NumPhases; ++k) {
				if (EDGES[j].hasphase[k]) {
					++numCon;
					variationReactiveConstraint1.add(IloRange(env, 0.0, IloInfinity));
					variationReactiveConstraint2.add(IloRange(env, 0.0, IloInfinity));
					for (int l = 0; l < NumPhases; ++l) {
						if (k != l && EDGES[j].hasphase[l]) {
							variationReactiveConstraint1[numCon].setLinearCoef(flowReactiveVariable[j*NumPhases + l], (PhaseVariation - 1.0) / (double) EDGES[j].NumPhases);
							variationReactiveConstraint2[numCon].setLinearCoef(flowReactiveVariable[j*NumPhases + l], (PhaseVariation + 1.0) / (double) EDGES[j].NumPhases);
						}
					}	
					variationReactiveConstraint1[numCon].setLinearCoef(flowReactiveVariable[j*NumPhases + k], ((double) EDGES[j].NumPhases + PhaseVariation - 1.0) / (double) EDGES[j].NumPhases);
					variationReactiveConstraint2[numCon].setLinearCoef(flowReactiveVariable[j*NumPhases + k], (PhaseVariation + 1.0 - (double) EDGES[j].NumPhases) / (double) EDGES[j].NumPhases);		
				}
			}
		}
	}
	// HARDEN CONSTRAINTS
	numCon = -1;
	cout << "HARDENED DISABLED: ";
	for (int j = 0; j < HARDENED_DISABLED.size(); ++j) {
		if (HARDENED_DISABLED[j].data && HardenedDisabled) {
			cout << HARDENED_DISABLED[j].id << " ";
			//lineHardenVariable[hashTableEdge[HARDENED_DISABLED[j].id]].setUB(0.0);
			++numCon;
			hardenConstraint.add(IloRange(env, -IloInfinity, 0.0));
			hardenConstraint[numCon].setLinearCoef(lineHardenVariable[hashTableEdge[HARDENED_DISABLED[j].id]], 1.0);
		}
	}
	cout << endl;
	// DAMAGE CONSTRAINTS
	if (ScenarioIndex != NumScenarios) { // BASE SCENARIO IS EVERYDAY = NO DAMAGE
		numCon = -1;
		cout << "DISABLED: ";
		for (int j = 0; j < DISABLED.size(); ++j) {
			//cout << "(" << DISABLED[j].id << "," << DISABLED[j].data << ") ";
			if (DISABLED[j].data) {
				cout << DISABLED[j].id << " ";
				// TO DAMAGE A LINE THAT ALREADY EXISTS, WE UPDATE ITS LOWER BOUND
				if (edgeExists[hashTableEdge[DISABLED[j].id]]) {
					lineUseVariable[hashTableEdge[DISABLED[j].id]].setLB(0);
				}
				++numCon;
				damageConstraint.add(IloRange(env, 0.0, 0.0));
				damageConstraint[numCon].setLinearCoef(lineUseVariable[hashTableEdge[DISABLED[j].id]], 1.0);
				if (supportable[hashTableEdge[DISABLED[j].id]]) {
					damageConstraint[numCon].setLinearCoef(lineHardenVariable[hashTableEdge[DISABLED[j].id]], -1.0);
				}
				//cout << "lineUseVariable " << hashTableEdge[DISABLED[j].id] << " disabled" << endl;
			}
		}
		cout << endl;
		NumDamagedEdges = numCon + 1;
	} else {
		NumDamagedEdges = 0;
	}
	// FOLLOWING CONSTRAINTS SHOULD BE UNITED UNDER A SINGLE LOOP
	// LOAD CONSTRAINTS
	numCon = -1;
	for (int j = 0; j < NODES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {	
			if (NODES[j].hasphase[k]) {
				++numCon;
				loadConstraint.add(IloRange(env, 0.0, 0.0));
				loadReactiveConstraint.add(IloRange(env, 0.0, 0.0));
				loadConstraint[numCon].setLinearCoef(loadVariable[j*NumPhases + k], -1.0);
				loadReactiveConstraint[numCon].setLinearCoef(loadReactiveVariable[j*NumPhases + k], -1.0);
				for (int l = 0; l < LOADS.size(); ++l) {
					if (hashTableVertex[LOADS[l].node_id] == j && LOADS[l].hasphase[k]) {
						loadConstraint[numCon].setLinearCoef(loadServeVariable[l], LOADS[l].realphase[k]);
						loadReactiveConstraint[numCon].setLinearCoef(loadServeVariable[l], LOADS[l].reactivephase[k]);
					}
				}	
			}
		}
	}
	// GENERATION CONSTRAINTS
	vector<vector<double> > MAXREALPHASE(NODES.size(), vector<double>(NumPhases, 0.0));
	vector<vector<double> > MAXREACTIVEPHASE(NODES.size(), vector<double>(NumPhases, 0.0));
	for (int j = 0; j < NODES.size(); ++j) {
		for (int k = 0; k < GENERATORS.size(); ++k) {
			for (int l = 0; l < NumPhases; ++l) {
				if (hashTableVertex[GENERATORS[k].node_id] == j && GENERATORS[k].hasphase[l]) {
					MAXREALPHASE[j][l] += GENERATORS[k].maxrealphase[l];
					MAXREACTIVEPHASE[j][l] += GENERATORS[k].maxreactivephase[l];
				}
			}
		}
	}
	vector<double> MAX_ADDED(GENERATORS.size(), 0.0);
	for (int j = 0; j < MAX_MICROGRID.size(); ++j) {
		MAX_ADDED[hashTableGenerator[MAX_MICROGRID[j].id]] = MAX_MICROGRID[j].data;
	}
	numCon = -1;
	for (int j = 0; j < NODES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			++numCon;
			generationConstraint.add(IloRange(env, -MAXREALPHASE[j][k], IloInfinity));
			generationConstraint[numCon].setLinearCoef(generatorVariable[j*NumPhases + k], -1.0);
			for (int l = 0; l < GENERATORS.size(); ++l) {	
				if (hashTableVertex[GENERATORS[l].node_id] == j && GENERATORS[l].hasphase[k]) {
					generationConstraint[numCon].setLinearCoef(facilityVariable[l], MAX_ADDED[l]);
				}
			}	
		}
	}
	numCon = -1;
	for (int j = 0; j < NODES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			++numCon;
			generationReactiveConstraint.add(IloRange(env, -MAXREACTIVEPHASE[j][k], IloInfinity));
			generationReactiveConstraint[numCon].setLinearCoef(generatorReactiveVariable[j*NumPhases + k], -1.0);
			for (int l = 0; l < GENERATORS.size(); ++l) {	
				if (hashTableVertex[GENERATORS[l].node_id] == j && GENERATORS[l].hasphase[k]) {
					generationReactiveConstraint[numCon].setLinearCoef(facilityVariable[l], MAX_ADDED[l]);
				}
			}	
		}
	}
	// BALANCE CONSTRAINTS		
	numCon = -1;
	for (int j = 0; j < NODES.size(); ++j) {
		//cout << NODES[j].hasphase[0] << " " << NODES[j].hasphase[1] << " " << NODES[j].hasphase[2] << " " << endl;
		for (int k = 0; k < NumPhases; ++k) {
			if (balancable[j][k]) {
				// NEED TO ADD FOR ALL NODES FOR FLOW CONSERVATION
				++numCon;
				balanceConstraint.add(IloRange(env, 0.0, 0.0));
				//list<int>::iterator ind = G[j].AdjList.begin();
				bool undeliverable = true;
				for (list<string>::iterator it = G[j].EdgeID.begin(); it != G[j].EdgeID.end(); ++it) {
					if (EDGES[hashTableEdge[*it]].hasphase[k]) {
						undeliverable = false;
						if (j == hashTableVertex[EDGES[hashTableEdge[*it]].node1id]) {//(j < (*ind)) { // Check which node edge contributes
							balanceConstraint[numCon].setLinearCoef(flowVariable[hashTableEdge[*it]*NumPhases + k], -1.0);
						} else {
							balanceConstraint[numCon].setLinearCoef(flowVariable[hashTableEdge[*it]*NumPhases + k], 1.0);
						}
					}
					//++ind;
				}
				if (undeliverable && NODES[j].hasphase[k]) { // UNDELIVERABLE NODE
					cout << "Node #" << j << " has undeliverable demand of " << NODES[j].demand[k] << " units on phase #" << k << endl;
				}
				if (NODES[j].hasphase[k]) {
					balanceConstraint[numCon].setLinearCoef(loadVariable[j*NumPhases + k], -1.0);
				}	
				balanceConstraint[numCon].setLinearCoef(generatorVariable[j * NumPhases + k], 1.0);
			}
		}
	}
	numCon = -1;
	for (int j = 0; j < NODES.size(); ++j) {
		//cout << NODES[j].hasphase[0] << " " << NODES[j].hasphase[1] << " " << NODES[j].hasphase[2] << " " << endl;
		for (int k = 0; k < NumPhases; ++k) {
			if (balancable[j][k]) {
				// NEED TO ADD FOR ALL NODES FOR FLOW CONSERVATION
				++numCon;
				balanceReactiveConstraint.add(IloRange(env, 0.0, 0.0));
				//list<int>::iterator ind = G[j].AdjList.begin();
				bool undeliverable = true;
				for (list<string>::iterator it = G[j].EdgeID.begin(); it != G[j].EdgeID.end(); ++it) {
					if (EDGES[hashTableEdge[*it]].hasphase[k]) {
						undeliverable = false;
						if (j == hashTableVertex[EDGES[hashTableEdge[*it]].node1id]) {//(j < (*ind)) { // Check which node edge contributes
							balanceReactiveConstraint[numCon].setLinearCoef(flowReactiveVariable[hashTableEdge[*it]*NumPhases + k], -1.0);
						} else {
							balanceReactiveConstraint[numCon].setLinearCoef(flowReactiveVariable[hashTableEdge[*it]*NumPhases + k], 1.0);
						}
					}
					//++ind;
				}
				if (undeliverable && NODES[j].hasphase[k]) { // UNDELIVERABLE NODE
					cout << "Node #" << j << " has undeliverable demand of " << NODES[j].demand[k] << " units on phase #" << k << endl;
				}
				if (NODES[j].hasphase[k]) {
					balanceReactiveConstraint[numCon].setLinearCoef(loadReactiveVariable[j*NumPhases + k], -1.0);
				}	
				balanceReactiveConstraint[numCon].setLinearCoef(generatorReactiveVariable[j * NumPhases + k], 1.0);
			}
		}
	}
	// MICROGRID CONSTRAINTS
	/*vector<double> MAX_ADDED(GENERATORS.size(), 0.0);
	for (int j = 0; j < MAX_MICROGRID.size(); ++j) {
		MAX_ADDED[hashTableGenerator[MAX_MICROGRID[j].id]] += MAX_MICROGRID[j].data;
	}
	numCon = -1;
	for (int j = 0; j < GENERATORS.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			if (GENERATORS[j].hasphase[k]) {
				++numCon;
				microgridConstraint.add(IloRange(env, -IloInfinity, 0.0));
				microgridConstraint[numCon].setLinearCoef(facilityVariable[j], -MAX_ADDED[j]);
				microgridConstraint[numCon].setLinearCoef(microgridVariable[j*NumPhases + k], 1.0);	
			}
		}
	}*/
	// CYCLE CONSTRAINTS
	// COPY EDGES ARE NEEDED TO BE ENUMERATED EN EACH CYCLE AS THEY MAY CARRY DISTINCT SET OF PHASES (they are not treated as cycles)
	numCon = -1;
	for (int j = 0; j < CYCLES.size(); ++j) {
		// FOLLOWING ENUMERATION IS EXPENSIVE AS THERE ARE 300k CYCLES OVER DISTINCT EDGES FOR SOME GRAPHS
		//vector<bool> include(EDGES.size(), false);
		//enumerateCycleConstraints(env, include, 0, cycleConstraint, lineUseVariable, switchUseVariable, CYCLES[j], CopyEdge, hashTableEdge, hashTableVertex, EDGES, numCon);

		cycleConstraint.add(IloRange(env, -IloInfinity, (double) CYCLES[j].size() - 1.0));
		for (int k = 0; k < EDGES.size(); ++k) {
			if (ismember(CYCLES[j], hashTableVertex[EDGES[k].node1id], hashTableVertex[EDGES[k].node2id])) {
				//cout << hashTableCopyEdge[EDGES[k].id] << endl;
				cycleConstraint[j].setLinearCoef(lineCycleVariable[hashTableCopyEdge[EDGES[k].id]], 1.0);
				cycleConstraint[j].setLinearCoef(switchCycleVariable[hashTableCopyEdge[EDGES[k].id]], -1.0);
			} 
		}
	}
	for (int j = 0; j < EDGES.size(); ++j) {
		cycleSubConstraint.add(IloRange(env, -IloInfinity, 0.0));
		cycleSubConstraint[j].setLinearCoef(lineUseVariable[j], 1.0);
		cycleSubConstraint[j].setLinearCoef(lineCycleVariable[hashTableCopyEdge[EDGES[j].id]], -1.0);
	}
	// SWITCHABLE CONSTRAINTS
	for (int j = 0; j < NumDistinctEdges; ++j) {
		switchableConstraint.add(IloRange(env, 0.0, IloInfinity));
		switchableConstraint[j].setLinearCoef(lineCycleVariable[j], 1.0);
		switchableConstraint[j].setLinearCoef(switchCycleVariable[j], -1.0);
	}
	// CRITICAL SERVE CONSTRAINT
	double CriticalDemand = 0.0;
	double CriticalReactiveDemand = 0.0;
	double TotalDemand = 0.0;
	double TotalReactiveDemand = 0.0;
	for (int j = 0; j < IS_CRITICAL_LOAD.size(); ++j) {
		if (IS_CRITICAL_LOAD[j].data) {
			CriticalDemand += vectorsum(LOADS[hashTableLoad[IS_CRITICAL_LOAD[j].id]].realphase);
			CriticalReactiveDemand += vectorsum(LOADS[hashTableLoad[IS_CRITICAL_LOAD[j].id]].reactivephase);
		}
	}
	for (int j = 0; j < LOADS.size(); ++j) {
		TotalDemand += vectorsum(LOADS[j].realphase);
		TotalReactiveDemand += vectorsum(LOADS[j].reactivephase);
	}
	criticalServeConstraint.add(IloRange(env, CriticalLoadMet * CriticalDemand, IloInfinity));
	criticalServeReactiveConstraint.add(IloRange(env, CriticalLoadMet * CriticalReactiveDemand, IloInfinity));
	loadServeConstraint.add(IloRange(env, LoadMet * TotalDemand, IloInfinity));
	loadServeReactiveConstraint.add(IloRange(env, LoadMet * TotalReactiveDemand, IloInfinity));

	//criticalServeConstraint.add(IloRange(env, 0.0, IloInfinity));
	//criticalServeReactiveConstraint.add(IloRange(env, 0.0, IloInfinity));
	//loadServeConstraint.add(IloRange(env, 0.0, IloInfinity));
	//loadServeReactiveConstraint.add(IloRange(env, 0.0, IloInfinity));

	//cout << "Critical Demand: " << CriticalDemand << " " << CriticalReactiveDemand << endl;
	//cout << "Total Demand: " << TotalDemand << " " << TotalReactiveDemand << endl;
	for (int j = 0; j < IS_CRITICAL_LOAD.size(); ++j) {
		if (IS_CRITICAL_LOAD[j].data) {
			for (int k = 0; k < NumPhases; ++k) {
				if (LOADS[hashTableLoad[IS_CRITICAL_LOAD[j].id]].hasphase[k]) {
					criticalServeConstraint[0].setLinearCoef(loadVariable[hashTableVertex[LOADS[hashTableLoad[IS_CRITICAL_LOAD[j].id]].node_id] * NumPhases + k], 1.0);
					criticalServeReactiveConstraint[0].setLinearCoef(loadReactiveVariable[hashTableVertex[LOADS[hashTableLoad[IS_CRITICAL_LOAD[j].id]].node_id] * NumPhases + k], 1.0);
				}
			}
		}
	}
	for (int j = 0; j < NODES.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			if (NODES[j].hasphase[k]) {
				loadServeConstraint[0].setLinearCoef(loadVariable[j * NumPhases + k], 1.0);
				loadServeReactiveConstraint[0].setLinearCoef(loadReactiveVariable[j * NumPhases + k], 1.0);
			}
		}
	}
	// SWITCH IMPLICATION CONSTRAINTS
	for (int j = 0; j < EDGES.size(); ++j) {
		switchImplicationConstraint1.add(IloRange(env, -1.0, IloInfinity));
		//switchImplicationConstraint2.add(IloRange(env, -IloInfinity, 3.0));
		switchImplicationConstraint1[j].setLinearCoef(switchUseVariable[j], 1.0);
		switchImplicationConstraint1[j].setLinearCoef(lineUseVariable[j], -1.0);
		switchImplicationConstraint1[j].setLinearCoef(switchCycleVariable[hashTableCopyEdge[EDGES[j].id]], -1.0);
		//switchImplicationConstraint2[j].setLinearCoef(switchUseVariable[j], 1.0);
		//switchImplicationConstraint2[j].setLinearCoef(lineUseVariable[j], 1.0);
		//switchImplicationConstraint2[j].setLinearCoef(switchCycleVariable[hashTableCopyEdge[EDGES[j].id]], 1.0);
	}


	GenerateObjectiveFunction(
		hashTableVertex,
		hashTableEdge,
		hashTableGenerator,
		hashTableLoad,
		NODES,
		EDGES,
		GENERATORS,
		LOADS,
		LINE_CONSTRUCTION_COST,
		MICROGRID_COST,
		MICROGRID_FIXED_COST,
		IS_CRITICAL_LOAD,
		MAX_MICROGRID,
		HARDEN_COST,
		LINE_SWITCH_COST,
		SubObjective,
		ObjectiveVector,
		lineUseVariable, // VARIABLES
		lineExistsVariable, // VARIABLES
		lineDirectionVariable,
		lineCycleVariable,
		lineHardenVariable,
		switchUseVariable,
		switchCycleVariable,
		flowVariable,
		flowReactiveVariable,
		voltageVariable,
		voltageOffsetVariable,
		loadVariable,
		loadReactiveVariable,
		loadServeVariable,
		generatorVariable,
		generatorReactiveVariable,
		//IloNumVarArray &microgridVariable,
		facilityVariable);



}

struct sortWithIndices {
	sortWithIndices() {ind = -1; val = 1e+10; tie = 0;}
	int ind;
	double val;
	int tie;
};

bool mycomperator(const sortWithIndices &a, const sortWithIndices &b) {
	return (fabs(a.val - b.val) < 1e-4 ? a.tie > b.tie : a.val < b.val);
}

bool ascendingcomperator(const double &a, const double &b) {
	return (a < b ? true : false);
}

template<class T>
T maxelement(const vector<T> &myvec) {
	T max = (T) 0.0;
	for (int i = 0; i < myvec.size(); ++i) {
		if (myvec[i] > max) {
			max = myvec[i];
		}
	}
	return max;
}

template<class T>
T minelement(const vector<T> &myvec) {
	T min = (T) 1e+10;
	for (int i = 0; i < myvec.size(); ++i) {
		if (myvec[i] < min) {
			min = myvec[i];
		}
	}
	return min;
}  


// SOLVES THE LP RELAXATION OF MASTER PROBLEM
void solveRelaxed(
	IloModel &MasterProblem,
	Solver &MasterSolver,
	IloNumVarArray &lineUseVariable, // VARIABLES
	IloNumVarArray &lineExistsVariable, // VARIABLES
	IloNumVarArray &lineDirectionVariable,
	IloNumVarArray &lineCycleVariable,
	IloNumVarArray &lineHardenVariable,
	IloNumVarArray &switchUseVariable,
	IloNumVarArray &switchCycleVariable,
	IloNumVarArray &flowVariable,
	IloNumVarArray &flowReactiveVariable,
	IloNumVarArray &voltageVariable,
	IloNumVarArray &voltageOffsetVariable,
	IloNumVarArray &loadVariable,
	IloNumVarArray &loadReactiveVariable,
	IloNumVarArray &loadServeVariable,
	IloNumVarArray &generatorVariable,
	IloNumVarArray &generatorReactiveVariable,
	//IloNumVarArray &microgridVariable,
	IloNumVarArray &facilityVariable,
	vector<IloNumVarArray> &clonelineUseVariable, // SUB PROBLEM VARIABLES
	vector<IloNumVarArray> &clonelineExistsVariable, // SUB PROBLEM VARIABLES
	vector<IloNumVarArray> &clonelineDirectionVariable,
	vector<IloNumVarArray> &clonelineCycleVariable,
	vector<IloNumVarArray> &clonelineHardenVariable,
	vector<IloNumVarArray> &cloneswitchUseVariable,
	vector<IloNumVarArray> &cloneswitchCycleVariable,
	vector<IloNumVarArray> &cloneflowVariable,
	vector<IloNumVarArray> &cloneflowReactiveVariable,
	vector<IloNumVarArray> &clonevoltageVariable,
	vector<IloNumVarArray> &clonevoltageOffsetVariable,
	vector<IloNumVarArray> &cloneloadVariable,
	vector<IloNumVarArray> &cloneloadReactiveVariable,
	vector<IloNumVarArray> &cloneloadServeVariable,
	vector<IloNumVarArray> &clonegeneratorVariable,
	vector<IloNumVarArray> &clonegeneratorReactiveVariable,
	//vector<IloNumVarArray> &clonemicrogridVariable,
	//vector<IloNumVarArray> &clonefacilityVariable,
	vector<double> &LPSolution,
	int NumScenarios,
	int NumEdges,
	int NumGenerators,
	int NumPhases) {

	IloEnv env = MasterSolver.cplex.getEnv();
	IloModel RelaxedProblem(env);
	RelaxedProblem.add(MasterProblem);

	RelaxedProblem.add(IloConversion(env, lineUseVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, lineExistsVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, lineDirectionVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, lineCycleVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, lineHardenVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, switchUseVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, switchCycleVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, flowVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, flowReactiveVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, voltageVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, voltageOffsetVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, loadVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, loadReactiveVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, loadServeVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, generatorVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, generatorReactiveVariable, ILOFLOAT));
	//RelaxedProblem.add(IloConversion(env, microgridVariable, ILOFLOAT));
	RelaxedProblem.add(IloConversion(env, facilityVariable, ILOFLOAT));

	for (int i = 0; i < clonelineUseVariable.size(); ++i) {
		RelaxedProblem.add(IloConversion(env, clonelineUseVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, clonelineExistsVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, clonelineDirectionVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, clonelineCycleVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, clonelineHardenVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, cloneswitchUseVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, cloneswitchCycleVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, cloneflowVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, cloneflowReactiveVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, clonevoltageVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, clonevoltageOffsetVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, cloneloadVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, cloneloadReactiveVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, cloneloadServeVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, clonegeneratorVariable[i], ILOFLOAT));
		RelaxedProblem.add(IloConversion(env, clonegeneratorReactiveVariable[i], ILOFLOAT));
		//RelaxedProblem.add(IloConversion(env, clonemicrogridVariable[i], ILOFLOAT));
		//RelaxedProblem.add(IloConversion(env, clonefacilityVariable[i], ILOFLOAT));
	}

	IloCplex LPsolver(RelaxedProblem);
	LPsolver.solve();

	IloNumArray lineusevals(env);
	IloNumArray lineswitchvals(env);
	IloNumArray linehardenvals(env);
	IloNumArray facilityvals(env);

	LPsolver.getValues(lineusevals, lineUseVariable);
	LPsolver.getValues(lineswitchvals, switchUseVariable);
	LPsolver.getValues(linehardenvals, lineHardenVariable);
	//LPsolver.getValues(microgridvals, microgridVariable);
	LPsolver.getValues(facilityvals, facilityVariable);

	for (int i = 0; i < NumEdges; ++i) {
		LPSolution[i] = lineusevals[i];
		LPSolution[i+NumEdges] = lineswitchvals[i];
		LPSolution[i+2*NumEdges] = linehardenvals[i];
	}
	for (int j = 0; j < NumGenerators; ++j) {
		LPSolution[3*NumEdges + j] = facilityvals[j];
	}


}

// THIS METHOD CAN BE USED TO GET OUT OF LOCAL MINIMA
bool VariableNeighborhoodDescent(
	double MIN_FEAS_OBJ,
	vector<double> &incumbentSolution,
	double &prevObjValue,
	vector<double> &LPSolution,
	vector<IloModel> &SubProblem,
	vector<Solver> &SubSolver,
	map<string, int> &hashTableVertex,
	map<string, int> &hashTableEdge,
	map<string, int> &hashTableGenerator,
	map<string, int> &hashTableLoad,
	vector<nodeData> &NODES,
	vector<edgeData> &EDGES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	map<int, lineCodeData> &LINECODES,
	vector<dataNode<double> > &LINE_CONSTRUCTION_COST,
	vector<dataNode<double> > &MICROGRID_COST,
	vector<dataNode<double> > &MICROGRID_FIXED_COST,
	vector<dataNode<bool> > &IS_CRITICAL_LOAD,
	vector<dataNode<double> > &MAX_MICROGRID,
	vector<dataNode<double> > &HARDEN_COST,
	vector<dataNode<double> > &LINE_SWITCH_COST,
	double LoadMet,
	double CriticalLoadMet,
	double PhaseVariation,
	vector<graph_vertex> &G,
	vector<vector<int> > &CYCLES,
	vector<IloNumVarArray> &lineUseVariable, // VARIABLES
	vector<IloNumVarArray> &lineExistsVariable, // VARIABLES
	vector<IloNumVarArray> &lineDirectionVariable,
	vector<IloNumVarArray> &lineCycleVariable,
	vector<IloNumVarArray> &lineHardenVariable,
	vector<IloNumVarArray> &switchUseVariable,
	vector<IloNumVarArray> &switchCycleVariable,
	vector<IloNumVarArray> &flowVariable,
	vector<IloNumVarArray> &flowReactiveVariable,
	vector<IloNumVarArray> &voltageVariable,
	vector<IloNumVarArray> &voltageOffsetVariable,
	vector<IloNumVarArray> &loadVariable,
	vector<IloNumVarArray> &loadReactiveVariable,
	vector<IloNumVarArray> &loadServeVariable,
	vector<IloNumVarArray> &generatorVariable,
	vector<IloNumVarArray> &generatorReactiveVariable,
	//vector<IloNumVarArray> &microgridVariable,
	vector<IloNumVarArray> &facilityVariable,
	const vector<int> &NumDamagedEdges,
	int NumPhases,
	int NumEdges,
	int NumGenerators,
	int NumScenarios,
	double discreteMicrogrid,
	string damage_str,
	string extent_str,
	string root_str,
	bool HardenedDisabled,
	double MAX_CPUTIME){

	IloEnv env = SubSolver[NumScenarios].cplex.getEnv();
	IloNumArray lineusevals(env);
	IloNumArray lineswitchvals(env);
	IloNumArray linehardenvals(env);
	IloNumArray facilityvals(env);

	bool retval = false;

	int ColSize = 3 * NumEdges + NumGenerators;

	// VARIABLE NEIGHBORHOOD DESCENT AROUND INCUMBENT
	double rhs = 1.0;
	double rhs_max = 5.0;//floor((double) ColSize / 100.0);
	time_t begin = clock();
	while (rhs < rhs_max) {
		//cout << "LOCAL NEIGHBORHOOD: " << rhs << endl;
		vector<int> isOne(ColSize, 0);
		for (int j = 0; j < ColSize; ++j) {
			if (fabs(incumbentSolution[j] - 1.0) < 1e-4) {
				isOne[j] = 1;
			}
		}
		IloRange DifferenceConstraint(env, -IloInfinity, (rhs - (double)vectorsum<int>(isOne))); 
		for (int j = 0; j < ColSize; ++j) {

			int fixIndex = j;
			int varIndex;
			if (fixIndex < NumEdges) {//lineUseVariable
				varIndex = fixIndex;
				if (isOne[fixIndex] == 1) {
					DifferenceConstraint.setLinearCoef(lineUseVariable[NumScenarios][varIndex], -1.0);
				} else {
					DifferenceConstraint.setLinearCoef(lineUseVariable[NumScenarios][varIndex], 1.0);
				}
			} else if (fixIndex < 2*NumEdges) {//switchUseVariable
				varIndex = fixIndex - NumEdges;
				if (isOne[fixIndex] == 1) {
					DifferenceConstraint.setLinearCoef(switchUseVariable[NumScenarios][varIndex], -1.0);
				} else {
					DifferenceConstraint.setLinearCoef(switchUseVariable[NumScenarios][varIndex], 1.0);
				}
			} else if (fixIndex < 3*NumEdges) {//lineHardenVariable
				varIndex = fixIndex - 2*NumEdges;
				if (isOne[fixIndex] == 1) {
					DifferenceConstraint.setLinearCoef(lineHardenVariable[NumScenarios][varIndex], -1.0);
				} else {
					DifferenceConstraint.setLinearCoef(lineHardenVariable[NumScenarios][varIndex], 1.0);
				}
			} else /*if (fixIndex < (3*NumEdges + NumPhases*NumGenerators)) {//microgridVariable
				
			} else*/ { //facilityVariable
				varIndex = fixIndex - 3*NumEdges - NumPhases*NumGenerators;
				if (isOne[fixIndex] == 1) {
					DifferenceConstraint.setLinearCoef(facilityVariable[NumScenarios][varIndex], -1.0);
				} else {
					DifferenceConstraint.setLinearCoef(facilityVariable[NumScenarios][varIndex], 1.0);
				}
			}
		}					
		SubProblem[NumScenarios].add(DifferenceConstraint);
		if ((MAX_CPUTIME - ((double) (clock() - begin) / (double) CLOCKS_PER_SEC)) < 60.0) {
			SubProblem[NumScenarios].remove(DifferenceConstraint);
			DifferenceConstraint.end();
			break;
		}
		SubSolver[NumScenarios].cplex.setParam(IloCplex::TiLim, MAX_CPUTIME - ((double) (clock() - begin) / (double) CLOCKS_PER_SEC));

		IloNumArray lineuseincumbent(env);
		IloNumArray lineswitchincumbent(env);
		IloNumArray linehardenincumbent(env);
		IloNumArray facilityincumbent(env);
		for (int j = 0; j < NumEdges; ++j) {
			lineuseincumbent.add(incumbentSolution[j]);
			lineswitchincumbent.add(incumbentSolution[j + NumEdges]);
			linehardenincumbent.add(incumbentSolution[j + 2*NumEdges]);
		}
		for (int j = 0; j < NumGenerators; ++j) {
			facilityincumbent.add(incumbentSolution[j + 3*NumEdges]);
		}

		SubSolver[NumScenarios].cplex.addMIPStart(lineUseVariable[NumScenarios], lineuseincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	
		SubSolver[NumScenarios].cplex.addMIPStart(switchUseVariable[NumScenarios], lineswitchincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	
		SubSolver[NumScenarios].cplex.addMIPStart(lineHardenVariable[NumScenarios], linehardenincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	
		//SubSolver[NumScenarios].cplex.addMIPStart(microgridVariable[NumScenarios], microgridincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	
		SubSolver[NumScenarios].cplex.addMIPStart(facilityVariable[NumScenarios], facilityincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	




		SubSolver[NumScenarios].cplex.solve();
		if (SubSolver[NumScenarios].cplex.getStatus() == IloAlgorithm::Feasible || SubSolver[NumScenarios].cplex.getStatus() == IloAlgorithm::Optimal) {
			if (SubSolver[NumScenarios].cplex.getObjValue() < prevObjValue - 1e-4) {
				retval = true;
				prevObjValue = SubSolver[NumScenarios].cplex.getObjValue();
				//cout << "OBJECTIVE: " << SubSolver[NumScenarios].cplex.getObjValue() << endl;
				SubSolver[NumScenarios].cplex.getValues(lineusevals, lineUseVariable[NumScenarios]);
				SubSolver[NumScenarios].cplex.getValues(lineswitchvals, switchUseVariable[NumScenarios]);
				SubSolver[NumScenarios].cplex.getValues(linehardenvals, lineHardenVariable[NumScenarios]);
				//SubSolver[NumScenarios].cplex.getValues(microgridvals, microgridVariable[NumScenarios]);
				SubSolver[NumScenarios].cplex.getValues(facilityvals, facilityVariable[NumScenarios]);

				// UPDATE INCUMBENT
				for (int i = 0; i < NumEdges; ++i) {
					incumbentSolution[i] = lineusevals[i];
					incumbentSolution[i+NumEdges] = lineswitchvals[i];
					incumbentSolution[i+2*NumEdges] = linehardenvals[i];
				}
				for (int j = 0; j < NumGenerators; ++j) {
					incumbentSolution[3*NumEdges + j] = facilityvals[j];
				}
				SubProblem[NumScenarios].remove(DifferenceConstraint);
				DifferenceConstraint.end();
				break;
			} else {
				rhs = rhs * 2.0;
			}
		} else {
			rhs = rhs * 2.0;
		}
		SubProblem[NumScenarios].remove(DifferenceConstraint);
		DifferenceConstraint.end();
	}
	return retval;

}



// FOLLOWING METHOD PERFORMS A VARIABLE NEIGHBORHOOD SEARCH OF MASTER PROBLEM
void VariableNeighborhoodSearchSolver(
	double MIN_FEAS_OBJ,
	vector<double> &incumbentSolution,
	vector<double> &LPSolution,
	vector<IloModel> &SubProblem,
	vector<Solver> &SubSolver,
	map<string, int> &hashTableVertex,
	map<string, int> &hashTableEdge,
	map<string, int> &hashTableGenerator,
	map<string, int> &hashTableLoad,
	vector<nodeData> &NODES,
	vector<edgeData> &EDGES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	map<int, lineCodeData> &LINECODES,
	vector<dataNode<double> > &LINE_CONSTRUCTION_COST,
	vector<dataNode<double> > &MICROGRID_COST,
	vector<dataNode<double> > &MICROGRID_FIXED_COST,
	vector<dataNode<bool> > &IS_CRITICAL_LOAD,
	vector<dataNode<double> > &MAX_MICROGRID,
	vector<dataNode<double> > &HARDEN_COST,
	vector<dataNode<double> > &LINE_SWITCH_COST,
	double LoadMet,
	double CriticalLoadMet,
	double PhaseVariation,
	vector<graph_vertex> &G,
	vector<vector<int> > &CYCLES,
	vector<IloNumVarArray> &lineUseVariable, // VARIABLES
	vector<IloNumVarArray> &lineExistsVariable, // VARIABLES
	vector<IloNumVarArray> &lineDirectionVariable,
	vector<IloNumVarArray> &lineCycleVariable,
	vector<IloNumVarArray> &lineHardenVariable,
	vector<IloNumVarArray> &switchUseVariable,
	vector<IloNumVarArray> &switchCycleVariable,
	vector<IloNumVarArray> &flowVariable,
	vector<IloNumVarArray> &flowReactiveVariable,
	vector<IloNumVarArray> &voltageVariable,
	vector<IloNumVarArray> &voltageOffsetVariable,
	vector<IloNumVarArray> &loadVariable,
	vector<IloNumVarArray> &loadReactiveVariable,
	vector<IloNumVarArray> &loadServeVariable,
	vector<IloNumVarArray> &generatorVariable,
	vector<IloNumVarArray> &generatorReactiveVariable,
	//vector<IloNumVarArray> &microgridVariable,
	vector<IloNumVarArray> &facilityVariable,
	vector<IloNumVarArray> &clonelineUseVariable,
	vector<IloNumVarArray> &clonelineExistsVariable,
	vector<IloNumVarArray> &clonelineDirectionVariable,
	vector<IloNumVarArray> &cloneswitchUseVariable,
	vector<IloNumVarArray> &clonegeneratorVariable,
	vector<IloNumVarArray> &clonegeneratorReactiveVariable,
	vector<IloNumVarArray> &clonevoltageVariable,
	vector<IloNumVarArray> &clonevoltageOffsetVariable,
	vector<IloNumVarArray> &cloneloadVariable,
	vector<IloNumVarArray> &cloneloadReactiveVariable,
	vector<IloNumVarArray> &cloneloadServeVariable,
	vector<IloNumArray> &cloneLineExistsVals,
	vector<IloNumArray> &cloneLineUseVals,
	vector<IloNumArray> &cloneLineDirectionVals,
	vector<IloNumArray> &cloneSwitchUseVals,
	vector<IloNumArray> &cloneGeneratorVals,
	vector<IloNumArray> &cloneGeneratorReactiveVals,
	vector<IloNumArray> &cloneLoadVals,
	vector<IloNumArray> &cloneLoadReactiveVals,
	vector<IloNumArray> &cloneLoadServeVals,
	vector<IloNumArray> &cloneVoltageVals,
	int bendersIter,
	const vector<int> &NumDamagedEdges,
	int NumPhases,
	int NumEdges,
	int NumGenerators,
	int NumScenarios,
	double discreteMicrogrid,
	string damage_str,
	string extent_str,
	string root_str,
	bool HardenedDisabled,
	double MAX_CPUTIME) {


	IloEnv env = SubSolver[NumScenarios].cplex.getEnv();

	int ColSize = 3 * NumEdges + NumGenerators;
	double d;
	double prevObjValue = MIN_FEAS_OBJ;
	bool reinitialize = false;
	int noimprovement = 0;
	int consecutiveReinitializationCounter = 0;
	int MAX_NO_IMPROVEMENT = 3;
	int MAX_REINITIALIZATION = 4;

	SubSolver[NumScenarios].cplex.setParam(IloCplex::ClockType, 1);
	SubSolver[NumScenarios].cplex.setParam(IloCplex::TiLim, 5.0 * 60.0);

	time_t begin = clock();
	while (((double) (clock() - begin) / (double) CLOCKS_PER_SEC) < MAX_CPUTIME) {

		if (consecutiveReinitializationCounter >= MAX_REINITIALIZATION) {
			if (false) {
			// FOLLOWING METHOD CAN BE USED TO ESCAPE LOCAL MINIMA
			/*if(VariableNeighborhoodDescent(
				MIN_FEAS_OBJ,
				incumbentSolution,
				prevObjValue,
				LPSolution,
				SubProblem,
				SubSolver,
				hashTableVertex,
				hashTableEdge,
				hashTableGenerator,
				hashTableLoad,
				NODES,
				EDGES,
				GENERATORS,
				LOADS,
				LINECODES,
				LINE_CONSTRUCTION_COST,
				MICROGRID_COST,
				MICROGRID_FIXED_COST,
				IS_CRITICAL_LOAD,
				MAX_MICROGRID,
				HARDEN_COST,
				LINE_SWITCH_COST,
				LoadMet,
				CriticalLoadMet,
				PhaseVariation,
				G,
				CYCLES,
				lineUseVariable, // VARIABLES
				lineExistsVariable, // VARIABLES
				lineDirectionVariable,
				lineCycleVariable,
				lineHardenVariable,
				switchUseVariable,
				switchCycleVariable,
				flowVariable,
				flowReactiveVariable,
				voltageVariable,
				voltageOffsetVariable,
				loadVariable,
				loadReactiveVariable,
				loadServeVariable,
				generatorVariable,
				generatorReactiveVariable,
				microgridVariable,
				facilityVariable,
				NumDamagedEdges,
				NumPhases,
				NumEdges,
				NumGenerators,
				NumScenarios,
				discreteMicrogrid,
				damage_str,
				extent_str,
				root_str,
				HardenedDisabled,
				5.0 * 60.0)) {
				consecutiveReinitializationCounter = 0;
				reinitialize = false;*/
			} else {
				break;
			}
		}

		noimprovement = 0;

		vector<sortWithIndices> delta(ColSize);
		for (int j = 0; j < ColSize; ++j) {
			delta[j].ind = j;
			delta[j].val = fabs(incumbentSolution[j] - LPSolution[j]);
		}

		if (reinitialize) {
			cout << "*** REINITIALIZATION ***" << endl;
			reinitialize = false;
			d = 0.5;
			std::random_shuffle(delta.begin(), delta.end());
			consecutiveReinitializationCounter++;
			for (int i = 0; i < 5; ++i) {
				cout << delta[i].ind << " ";
			}
			cout << endl;
		} else {
			d = 2.0;
			std::sort(delta.begin(), delta.end(), mycomperator);
		}		

		int NumNonZero = 0;
		for (int j = 0; j < ColSize; ++j) {
			if (delta[j].val > 1e-4) {
				NumNonZero++;
			}
		}
		int k_step = int((double)NumNonZero / d);
		int k_neighborhood = ColSize - k_step;
		
		do {


			vector<vector<double> > fixedBounds(0);
			for (int j = 0; j < k_neighborhood; ++j) {

				fixedBounds.push_back(vector<double>(2));
				int fixIndex = delta[j].ind;
				int varIndex;
				if (fixIndex < NumEdges) {//lineUseVariable
					varIndex = fixIndex;
					fixedBounds[j][0] = lineUseVariable[NumScenarios][varIndex].getLB();
					fixedBounds[j][1] = lineUseVariable[NumScenarios][varIndex].getUB();
					lineUseVariable[NumScenarios][varIndex].setLB(round(incumbentSolution[fixIndex]));
					lineUseVariable[NumScenarios][varIndex].setUB(round(incumbentSolution[fixIndex]));
				} else if (fixIndex < 2*NumEdges) {//switchUseVariable
					varIndex = fixIndex - NumEdges;
					fixedBounds[j][0] = switchUseVariable[NumScenarios][varIndex].getLB();
					fixedBounds[j][1] = switchUseVariable[NumScenarios][varIndex].getUB();
					switchUseVariable[NumScenarios][varIndex].setLB(round(incumbentSolution[fixIndex]));
					switchUseVariable[NumScenarios][varIndex].setUB(round(incumbentSolution[fixIndex]));
				} else if (fixIndex < 3*NumEdges) {//lineHardenVariable
					varIndex = fixIndex - 2*NumEdges;
					fixedBounds[j][0] = lineHardenVariable[NumScenarios][varIndex].getLB();
					fixedBounds[j][1] = lineHardenVariable[NumScenarios][varIndex].getUB();
					lineHardenVariable[NumScenarios][varIndex].setLB(round(incumbentSolution[fixIndex]));
					lineHardenVariable[NumScenarios][varIndex].setUB(round(incumbentSolution[fixIndex]));
				} else /*if (fixIndex < (3*NumEdges + NumPhases*NumGenerators)) {//microgridVariable
					varIndex = fixIndex - 3*NumEdges;
					fixedBounds[j][0] = microgridVariable[NumScenarios][varIndex].getLB();
					fixedBounds[j][1] = microgridVariable[NumScenarios][varIndex].getUB();
					microgridVariable[NumScenarios][varIndex].setLB(incumbentSolution[fixIndex] > 1e-4 ? incumbentSolution[fixIndex] : 0.0);
					microgridVariable[NumScenarios][varIndex].setUB(incumbentSolution[fixIndex] > 1e-4 ? incumbentSolution[fixIndex] : 0.0);
				} else*/ { //facilityVariable
					varIndex = fixIndex - 3*NumEdges;
					fixedBounds[j][0] = facilityVariable[NumScenarios][varIndex].getLB();
					fixedBounds[j][1] = facilityVariable[NumScenarios][varIndex].getUB();
					facilityVariable[NumScenarios][varIndex].setLB(round(incumbentSolution[fixIndex]));
					facilityVariable[NumScenarios][varIndex].setUB(round(incumbentSolution[fixIndex]));
				} 
				
			}	
			
			IloNumArray lineuseincumbent(env);
			IloNumArray lineswitchincumbent(env);
			IloNumArray linehardenincumbent(env);
			IloNumArray facilityincumbent(env);
			for (int j = 0; j < NumEdges; ++j) {
				lineuseincumbent.add(incumbentSolution[j]);
				lineswitchincumbent.add(incumbentSolution[j + NumEdges]);
				linehardenincumbent.add(incumbentSolution[j + 2*NumEdges]);
			}
			for (int j = 0; j < NumGenerators; ++j) {
				facilityincumbent.add(incumbentSolution[j + 3*NumEdges]);
			}

			SubSolver[NumScenarios].cplex.addMIPStart(lineUseVariable[NumScenarios], lineuseincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	
			SubSolver[NumScenarios].cplex.addMIPStart(switchUseVariable[NumScenarios], lineswitchincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	
			SubSolver[NumScenarios].cplex.addMIPStart(lineHardenVariable[NumScenarios], linehardenincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	
			//SubSolver[NumScenarios].cplex.addMIPStart(microgridVariable[NumScenarios], microgridincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	
			SubSolver[NumScenarios].cplex.addMIPStart(facilityVariable[NumScenarios], facilityincumbent, IloCplex::MIPStartAuto, "secondMIPStart");	


			SubSolver[NumScenarios].cplex.solve();


			cout << "NEIGHBORHOOD SIZE: " << k_neighborhood << " FREED NEIGHBORHOOD: " << k_step << " STATUS: " << SubSolver[NumScenarios].cplex.getStatus() << endl;
			if ((SubSolver[NumScenarios].cplex.getStatus() == IloAlgorithm::Optimal || SubSolver[NumScenarios].cplex.getStatus() == IloAlgorithm::Feasible)) {
				// COMPUTE ALL TREES IN THE CURRENT GLOBAL SOLUTION:
				cloneLineExistsVals.clear();
				cloneLineUseVals.clear();
				cloneLineDirectionVals.clear();
				cloneSwitchUseVals.clear();
				cloneGeneratorVals.clear();
				cloneGeneratorReactiveVals.clear();
				cloneLoadVals.clear();
				cloneLoadReactiveVals.clear();
				cloneLoadServeVals.clear();
				cloneVoltageVals.clear();
				for (int j = 0; j < (bendersIter+1); ++j) {
					cloneLineExistsVals.push_back(IloNumArray(env));
					cloneLineUseVals.push_back(IloNumArray(env));
					cloneLineDirectionVals.push_back(IloNumArray(env));
					cloneSwitchUseVals.push_back(IloNumArray(env));
					cloneGeneratorVals.push_back(IloNumArray(env));
					cloneGeneratorReactiveVals.push_back(IloNumArray(env));
					cloneLoadVals.push_back(IloNumArray(env));
					cloneLoadReactiveVals.push_back(IloNumArray(env));
					cloneLoadServeVals.push_back(IloNumArray(env));
					cloneVoltageVals.push_back(IloNumArray(env));
					//cout << clonelineUseVariable[j].getSize() << endl;
					SubSolver[NumScenarios].cplex.getValues(cloneLineExistsVals[j], clonelineExistsVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLineUseVals[j], clonelineUseVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLineDirectionVals[j], clonelineDirectionVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneSwitchUseVals[j], cloneswitchUseVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneGeneratorVals[j], clonegeneratorVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneGeneratorReactiveVals[j], clonegeneratorReactiveVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLoadVals[j], cloneloadVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLoadReactiveVals[j], cloneloadReactiveVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLoadServeVals[j], cloneloadServeVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneVoltageVals[j], clonevoltageVariable[j]);
				}
			}
			

			if ((SubSolver[NumScenarios].cplex.getStatus() == IloAlgorithm::Optimal || SubSolver[NumScenarios].cplex.getStatus() == IloAlgorithm::Feasible) && SubSolver[NumScenarios].cplex.getObjValue() < prevObjValue - 1e-4) {


				consecutiveReinitializationCounter = 0;
				prevObjValue = SubSolver[NumScenarios].cplex.getObjValue();
				env = SubSolver[NumScenarios].cplex.getEnv();
				IloNumArray lineusevals(env);
				IloNumArray lineswitchvals(env);
				IloNumArray linehardenvals(env);
				IloNumArray facilityvals(env);

				SubSolver[NumScenarios].cplex.getValues(lineusevals, lineUseVariable[NumScenarios]);
				SubSolver[NumScenarios].cplex.getValues(lineswitchvals, switchUseVariable[NumScenarios]);
				SubSolver[NumScenarios].cplex.getValues(linehardenvals, lineHardenVariable[NumScenarios]);
				//SubSolver[NumScenarios].cplex.getValues(microgridvals, microgridVariable[NumScenarios]);
				SubSolver[NumScenarios].cplex.getValues(facilityvals, facilityVariable[NumScenarios]);




				cout << "CURRENT STATUS: " << SubSolver[NumScenarios].cplex.getStatus() << endl;
				cout << "CURRENT OBJECTIVE VALUE: " << SubSolver[NumScenarios].cplex.getObjValue() << endl;
				cout << "CPUTIME: " << (double) (clock() - begin) / (double) CLOCKS_PER_SEC <<  endl;

				int numEdgesUsed = 0;
				int numSwitchesUsed = 0;
				int numHardenUsed = 0;
				int numFacilityUsed = 0;
				for (int j = 0; j < NumEdges; ++j) {
					if (fabs((double)lineusevals[j] - 1.0) < 1e-4) numEdgesUsed++;
					if (fabs((double)lineswitchvals[j] - 1.0) < 1e-4) numSwitchesUsed++;
					if (fabs((double)linehardenvals[j] - 1.0) < 1e-4) numHardenUsed++;
				}
				for (int j = 0; j < NumGenerators; ++j) {
					if (fabs((double)facilityvals[j] - 1.0) < 1e-4) numFacilityUsed++;
				}
				cout << "NUM EDGES: " << numEdgesUsed << " NUM SWITCHES: " << numSwitchesUsed << " NUM HARDEN: " << numHardenUsed << " NUM FACILITY: " << numFacilityUsed << endl;


				for (int j = 0; j < k_neighborhood; ++j) {

					int fixIndex = delta[j].ind;
					int varIndex;
					if (fixIndex < NumEdges) {//lineUseVariable
						varIndex = fixIndex;
						lineUseVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						lineUseVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} else if (fixIndex < 2*NumEdges) {//switchUseVariable
						varIndex = fixIndex - NumEdges;
						switchUseVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						switchUseVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} else if (fixIndex < 3*NumEdges) {//lineHardenVariable
						varIndex = fixIndex - 2*NumEdges;
						lineHardenVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						lineHardenVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} else /*if (fixIndex < (3*NumEdges + NumPhases*NumGenerators)) {//microgridVariable
						varIndex = fixIndex - 3*NumEdges;
						microgridVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						microgridVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} else*/ { //facilityVariable
						varIndex = fixIndex - 3*NumEdges;
						facilityVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						facilityVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} 
					
				}
				// UPDATE INCUMBENT
				for (int i = 0; i < NumEdges; ++i) {
					incumbentSolution[i] = lineusevals[i];
					incumbentSolution[i+NumEdges] = lineswitchvals[i];
					incumbentSolution[i+2*NumEdges] = linehardenvals[i];
				}
				for (int j = 0; j < NumGenerators; ++j) {
					incumbentSolution[3*NumEdges + j] = facilityvals[j];
				}		
	
				break;

			} else {
				for (int j = 0; j < k_neighborhood; ++j) {

					int fixIndex = delta[j].ind;
					int varIndex;
					if (fixIndex < NumEdges) {//lineUseVariable
						varIndex = fixIndex;
						lineUseVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						lineUseVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} else if (fixIndex < 2*NumEdges) {//switchUseVariable
						varIndex = fixIndex - NumEdges;
						switchUseVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						switchUseVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} else if (fixIndex < 3*NumEdges) {//lineHardenVariable
						varIndex = fixIndex - 2*NumEdges;
						lineHardenVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						lineHardenVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} else /*if (fixIndex < (3*NumEdges + NumPhases*NumGenerators)) {//microgridVariable
						varIndex = fixIndex - 3*NumEdges;
						microgridVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						microgridVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} else*/ { //facilityVariable
						varIndex = fixIndex - 3*NumEdges;
						facilityVariable[NumScenarios][varIndex].setLB(fixedBounds[j][0]);
						facilityVariable[NumScenarios][varIndex].setUB(fixedBounds[j][1]);
					} 
					
				}

				++noimprovement;
				if (noimprovement > MAX_NO_IMPROVEMENT) {
					reinitialize = true;
					break;
				}
				if (true) {// (k_neighborhood - k_step > ColSize - NumNonZero) {
					k_step = (int) (floor((double) k_step / 2.0) > 1.0 ? floor((double) k_step / 2.0) : 1.0);
				}
				k_neighborhood = k_neighborhood - k_step;
			}

		} while (true);
	
	}


	time_t end = clock();

	cout << "STATUS: Feasible OBJECTIVE: " << prevObjValue << endl;
	cout << "CPUTIME: " << (double) (end - begin) / (double) CLOCKS_PER_SEC <<  endl;

}

// FOLLOWING METHOD CAN BE USED TO SOLVE THE FULL MODEL WITHOUT DECOMPOSITION
void VariableNeighborhoodSearchMIP(
	vector<double> &ObjectiveVector,
	double MIN_FEAS_OBJ,
	vector<double> &MIN_FEAS_SOL,
	vector<IloModel> &SubProblem,
	vector<Solver> &SubSolver,
	map<string, int> &hashTableVertex,
	map<string, int> &hashTableEdge,
	map<string, int> &hashTableGenerator,
	map<string, int> &hashTableLoad,
	vector<nodeData> &NODES,
	vector<edgeData> &EDGES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	map<int, lineCodeData> &LINECODES,
	vector<dataNode<double> > &LINE_CONSTRUCTION_COST,
	vector<dataNode<double> > &MICROGRID_COST,
	vector<dataNode<double> > &MICROGRID_FIXED_COST,
	vector<dataNode<bool> > &IS_CRITICAL_LOAD,
	vector<dataNode<double> > &MAX_MICROGRID,
	vector<dataNode<double> > &HARDEN_COST,
	vector<dataNode<double> > &LINE_SWITCH_COST,
	double LoadMet,
	double CriticalLoadMet,
	double PhaseVariation,
	vector<graph_vertex> &G,
	vector<vector<int> > &CYCLES,
	vector<IloObjective> &SubObjective, // OBJECTIVE 
	vector<IloRangeArray> &lineExistsConstraint, // CONSTRAINTS
	vector<IloRangeArray> &lineConstraint1, // CONSTRAINTS
	vector<IloRangeArray> &lineConstraint2,
	vector<IloRangeArray> &lineReactiveConstraint1,
	vector<IloRangeArray> &lineReactiveConstraint2,
	vector<IloRangeArray> &directionConstraint,
	vector<IloRangeArray> &switchConstraint1,
	vector<IloRangeArray> &switchConstraint2,
	vector<IloRangeArray> &switchReactiveConstraint1,
	vector<IloRangeArray> &switchReactiveConstraint2,
	vector<IloRangeArray> &variationConstraint1,
	vector<IloRangeArray> &variationConstraint2,
	vector<IloRangeArray> &variationReactiveConstraint1,
	vector<IloRangeArray> &variationReactiveConstraint2,
	vector<IloRangeArray> &damageConstraint,
	vector<IloRangeArray> &loadConstraint,
	vector<IloRangeArray> &loadReactiveConstraint,
	vector<IloRangeArray> &generationConstraint,
	vector<IloRangeArray> &generationReactiveConstraint,
	vector<IloRangeArray> &balanceConstraint,
	vector<IloRangeArray> &balanceReactiveConstraint,
	vector<IloRangeArray> &microgridConstraint,
	vector<IloRangeArray> &cycleConstraint,
	vector<IloRangeArray> &cycleSubConstraint,
	vector<IloRangeArray> &switchableConstraint,
	vector<IloRangeArray> &criticalServeConstraint,
	vector<IloRangeArray> &criticalServeReactiveConstraint,
	vector<IloRangeArray> &loadServeConstraint,
	vector<IloRangeArray> &loadServeReactiveConstraint,
	vector<IloRangeArray> &switchImplicationConstraint1,
	//vector<IloRangeArray> &switchImplicationConstraint2,
	vector<IloRangeArray> &hardenConstraint,
	vector<IloRangeArray> &voltageConstraint,
	vector<IloRangeArray> &voltageOffsetConstraint,
	vector<IloRangeArray> &voltageVariationConstraint1,
	vector<IloRangeArray> &voltageVariationConstraint2,
	vector<IloNumVarArray> &lineUseVariable, // VARIABLES
	vector<IloNumVarArray> &lineExistsVariable, // VARIABLES
	vector<IloNumVarArray> &lineDirectionVariable,
	vector<IloNumVarArray> &lineCycleVariable,
	vector<IloNumVarArray> &lineHardenVariable,
	vector<IloNumVarArray> &switchUseVariable,
	vector<IloNumVarArray> &switchCycleVariable,
	vector<IloNumVarArray> &flowVariable,
	vector<IloNumVarArray> &flowReactiveVariable,
	vector<IloNumVarArray> &voltageVariable,
	vector<IloNumVarArray> &voltageOffsetVariable,
	vector<IloNumVarArray> &loadVariable,
	vector<IloNumVarArray> &loadReactiveVariable,
	vector<IloNumVarArray> &loadServeVariable,
	vector<IloNumVarArray> &generatorVariable,
	vector<IloNumVarArray> &generatorReactiveVariable,
	//vector<IloNumVarArray> &microgridVariable,
	vector<IloNumVarArray> &facilityVariable,
	const vector<int> &NumDamagedEdges,
	int NumPhases,
	int NumEdges,
	int NumGenerators,
	int NumScenarios,
	double discreteMicrogrid,
	string damage_str,
	string extent_str,
	string root_str,
	bool HardenedDisabled,
	vector<vector<dataNode<bool> > > &DISABLED,
	vector<vector<dataNode<bool> > > &HARDENED_DISABLED,
	bool CHECK_RESULTS,
	double newImpedanceMultiplier,
	bool LDFIndicator) {
	
	bool isInfeasible = true;	
	int numCon = -1;	

	int ColSize = 3 * NumEdges + NumGenerators;

	IloEnv env = SubSolver[NumScenarios].cplex.getEnv();

	IloRangeArray InfeasibilityCuts(env);
	IloNumVarArray InfeasibilityVars(env);
	int numInfeasibilityCuts = -1;

	IloRangeArray BaseLineVarConstraint(env);
	IloRangeArray BaseSwitchVarConstraint(env);
	IloRangeArray BaseHardenVarConstraint(env);
	IloRangeArray BaseMicrogridVarConstraint(env);
	IloRangeArray BaseFacilityVarConstraint(env);

	SubProblem[NumScenarios].add(BaseLineVarConstraint);
	SubProblem[NumScenarios].add(BaseSwitchVarConstraint);
	SubProblem[NumScenarios].add(BaseHardenVarConstraint);
	SubProblem[NumScenarios].add(BaseMicrogridVarConstraint);
	SubProblem[NumScenarios].add(BaseFacilityVarConstraint);

	IloObjective cloneSubObjective = IloObjective(env);

	vector<IloRangeArray> clonelineExistsConstraint(0);
	vector<IloRangeArray> clonelineConstraint1(0);
	vector<IloRangeArray> clonelineConstraint2(0);
	vector<IloRangeArray> clonelineReactiveConstraint1(0);
	vector<IloRangeArray> clonelineReactiveConstraint2(0);
	vector<IloRangeArray> clonedirectionConstraint(0);
	vector<IloRangeArray> cloneswitchConstraint1(0);
	vector<IloRangeArray> cloneswitchConstraint2(0);
	vector<IloRangeArray> cloneswitchReactiveConstraint1(0);
	vector<IloRangeArray> cloneswitchReactiveConstraint2(0);
	vector<IloRangeArray> clonevariationConstraint1(0);
	vector<IloRangeArray> clonevariationConstraint2(0);
	vector<IloRangeArray> clonevariationReactiveConstraint1(0);
	vector<IloRangeArray> clonevariationReactiveConstraint2(0);
	vector<IloRangeArray> clonedamageConstraint(0);
	vector<IloRangeArray> cloneloadConstraint(0);
	vector<IloRangeArray> cloneloadReactiveConstraint(0);
	vector<IloRangeArray> clonegenerationConstraint(0);
	vector<IloRangeArray> clonegenerationReactiveConstraint(0);
	vector<IloRangeArray> clonebalanceConstraint(0);
	vector<IloRangeArray> clonebalanceReactiveConstraint(0);
	vector<IloRangeArray> clonemicrogridConstraint(0);
	vector<IloRangeArray> clonecycleConstraint(0);
	vector<IloRangeArray> clonecycleSubConstraint(0);
	vector<IloRangeArray> cloneswitchableConstraint(0);
	vector<IloRangeArray> clonecriticalServeConstraint(0);
	vector<IloRangeArray> clonecriticalServeReactiveConstraint(0);
	vector<IloRangeArray> cloneloadServeConstraint(0);
	vector<IloRangeArray> cloneloadServeReactiveConstraint(0);
	vector<IloRangeArray> cloneswitchImplicationConstraint1(0);
	//vector<IloRangeArray> cloneswitchImplicationConstraint2(0);
	vector<IloRangeArray> clonehardenConstraint(0);
	vector<IloRangeArray> clonevoltageConstraint(0);
	vector<IloRangeArray> clonevoltageOffsetConstraint(0);
	vector<IloRangeArray> clonevoltageVariationConstraint1(0);
	vector<IloRangeArray> clonevoltageVariationConstraint2(0);

	vector<IloNumVarArray> clonelineUseVariable(0);
	vector<IloNumVarArray> clonelineExistsVariable(0);
	vector<IloNumVarArray> clonelineDirectionVariable(0);
	vector<IloNumVarArray> clonelineCycleVariable(0);
	vector<IloNumVarArray> clonelineHardenVariable(0);
	vector<IloNumVarArray> cloneswitchUseVariable(0);
	vector<IloNumVarArray> cloneswitchCycleVariable(0);
	vector<IloNumVarArray> cloneflowVariable(0);
	vector<IloNumVarArray> cloneflowReactiveVariable(0);
	vector<IloNumVarArray> clonevoltageVariable(0);
	vector<IloNumVarArray> clonevoltageOffsetVariable(0);
	vector<IloNumVarArray> cloneloadVariable(0);
	vector<IloNumVarArray> cloneloadReactiveVariable(0);
	vector<IloNumVarArray> cloneloadServeVariable(0);
	vector<IloNumVarArray> clonegeneratorVariable(0);
	vector<IloNumVarArray> clonegeneratorReactiveVariable(0);
	//vector<IloNumVarArray> clonemicrogridVariable(0);
	vector<IloNumVarArray> clonefacilityVariable(0);

	//vector<double> basesolution(ColSize, 0.0);
	//vector<vector<double> > subsolution(NumScenarios, vector<double>(ColSize, 0.0));

	//ofstream solutionFile("basesolutions.txt");

	cout << "NUM DAMAGED STATS:" << endl;
	cout << "Max element = " << maxelement<int>(NumDamagedEdges) << endl;
	cout << "Min element = " << minelement<int>(NumDamagedEdges) << endl;
	/*for (int i = 0; i < NumDamagedEdges.size(); ++i) {
		cout << NumDamagedEdges[i] << " ";
	}
	cout << endl;*/

	vector<vector<double> > MAXREALPHASE(NODES.size(), vector<double>(3, 0.0));
	for (int j = 0; j < NODES.size(); ++j) {
		for (int k = 0; k < GENERATORS.size(); ++k) {
			for (int l = 0; l < NumPhases; ++l) {
				if (hashTableVertex[GENERATORS[k].node_id] == j && GENERATORS[k].hasphase[l]) {
					MAXREALPHASE[j][l] += GENERATORS[k].maxrealphase[l];
				}
			}
		}
	}
	vector<double> MAX_ADDED(GENERATORS.size(), 0.0);
	for (int j = 0; j < MAX_MICROGRID.size(); ++j) {
		MAX_ADDED[hashTableGenerator[MAX_MICROGRID[j].id]] = MAX_MICROGRID[j].data;
	}
	int bendersIter = -1;
	for (int i = 0; i < NumScenarios; ++i) {
		++bendersIter;
		
		env = SubSolver[NumScenarios].cplex.getEnv();
		int ScenarioIndex = i;

		int numDamagedEdges;

		clonelineExistsConstraint.push_back(IloRangeArray(env));
		clonelineConstraint1.push_back(IloRangeArray(env));
		clonelineConstraint2.push_back(IloRangeArray(env));
		clonelineReactiveConstraint1.push_back(IloRangeArray(env));
		clonelineReactiveConstraint2.push_back(IloRangeArray(env));
		clonedirectionConstraint.push_back(IloRangeArray(env));
		cloneswitchConstraint1.push_back(IloRangeArray(env));
		cloneswitchConstraint2.push_back(IloRangeArray(env));
		cloneswitchReactiveConstraint1.push_back(IloRangeArray(env));
		cloneswitchReactiveConstraint2.push_back(IloRangeArray(env));
		clonevariationConstraint1.push_back(IloRangeArray(env));
		clonevariationConstraint2.push_back(IloRangeArray(env));
		clonevariationReactiveConstraint1.push_back(IloRangeArray(env));
		clonevariationReactiveConstraint2.push_back(IloRangeArray(env));
		clonedamageConstraint.push_back(IloRangeArray(env));
		cloneloadConstraint.push_back(IloRangeArray(env));
		cloneloadReactiveConstraint.push_back(IloRangeArray(env));
		clonegenerationConstraint.push_back(IloRangeArray(env));
		clonegenerationReactiveConstraint.push_back(IloRangeArray(env));
		clonebalanceConstraint.push_back(IloRangeArray(env));
		clonebalanceReactiveConstraint.push_back(IloRangeArray(env));
		clonemicrogridConstraint.push_back(IloRangeArray(env));
		clonecycleConstraint.push_back(IloRangeArray(env));
		clonecycleSubConstraint.push_back(IloRangeArray(env));
		cloneswitchableConstraint.push_back(IloRangeArray(env));
		clonecriticalServeConstraint.push_back(IloRangeArray(env));
		clonecriticalServeReactiveConstraint.push_back(IloRangeArray(env));
		cloneloadServeConstraint.push_back(IloRangeArray(env));
		cloneloadServeReactiveConstraint.push_back(IloRangeArray(env));
		cloneswitchImplicationConstraint1.push_back(IloRangeArray(env));
		//cloneswitchImplicationConstraint2.push_back(IloRangeArray(env));
		clonehardenConstraint.push_back(IloRangeArray(env));
		clonevoltageConstraint.push_back(IloRangeArray(env));
		clonevoltageOffsetConstraint.push_back(IloRangeArray(env));
		clonevoltageVariationConstraint1.push_back(IloRangeArray(env));
		clonevoltageVariationConstraint2.push_back(IloRangeArray(env));

		clonelineUseVariable.push_back(IloNumVarArray(env));
		clonelineExistsVariable.push_back(IloNumVarArray(env));
		clonelineDirectionVariable.push_back(IloNumVarArray(env));
		clonelineCycleVariable.push_back(IloNumVarArray(env));
		clonelineHardenVariable.push_back(IloNumVarArray(env));
		cloneswitchUseVariable.push_back(IloNumVarArray(env));
		cloneswitchCycleVariable.push_back(IloNumVarArray(env));
		cloneflowVariable.push_back(IloNumVarArray(env));
		cloneflowReactiveVariable.push_back(IloNumVarArray(env));
		clonevoltageVariable.push_back(IloNumVarArray(env));
		clonevoltageOffsetVariable.push_back(IloNumVarArray(env));
		cloneloadVariable.push_back(IloNumVarArray(env));
		cloneloadReactiveVariable.push_back(IloNumVarArray(env));
		cloneloadServeVariable.push_back(IloNumVarArray(env));
		clonegeneratorVariable.push_back(IloNumVarArray(env));
		clonegeneratorReactiveVariable.push_back(IloNumVarArray(env));
		//clonemicrogridVariable.push_back(IloNumVarArray(env));
		clonefacilityVariable.push_back(IloNumVarArray(env));


		vector<double> redObjectiveVector(ColSize,0.0);

		PopulateScenarioData(
			SubSolver[NumScenarios],
			hashTableVertex,
			hashTableEdge,
			hashTableGenerator,
			hashTableLoad,
			NODES,
			EDGES,
			GENERATORS,
			LOADS,
			LINECODES,
			LINE_CONSTRUCTION_COST,
			MICROGRID_COST,
			MICROGRID_FIXED_COST,
			IS_CRITICAL_LOAD,
			MAX_MICROGRID,
			HARDEN_COST,
			LINE_SWITCH_COST,
			LoadMet,
			CriticalLoadMet,
			PhaseVariation,
			G,
			CYCLES,
			cloneSubObjective,
			redObjectiveVector,
			clonelineExistsConstraint[bendersIter], // CONSTRAINTS
			clonelineConstraint1[bendersIter], // CONSTRAINTS
			clonelineConstraint2[bendersIter],
			clonelineReactiveConstraint1[bendersIter], // CONSTRAINTS
			clonelineReactiveConstraint2[bendersIter],
			clonedirectionConstraint[bendersIter],
			cloneswitchConstraint1[bendersIter],
			cloneswitchConstraint2[bendersIter],
			cloneswitchReactiveConstraint1[bendersIter],
			cloneswitchReactiveConstraint2[bendersIter],
			clonevariationConstraint1[bendersIter],
			clonevariationConstraint2[bendersIter],
			clonevariationReactiveConstraint1[bendersIter],
			clonevariationReactiveConstraint2[bendersIter],
			clonedamageConstraint[bendersIter],
			cloneloadConstraint[bendersIter],
			cloneloadReactiveConstraint[bendersIter],
			clonegenerationConstraint[bendersIter],
			clonegenerationReactiveConstraint[bendersIter],
			clonebalanceConstraint[bendersIter],
			clonebalanceReactiveConstraint[bendersIter],
			clonemicrogridConstraint[bendersIter],
			clonecycleConstraint[bendersIter],
			clonecycleSubConstraint[bendersIter],
			cloneswitchableConstraint[bendersIter],
			clonecriticalServeConstraint[bendersIter],
			clonecriticalServeReactiveConstraint[bendersIter],
			cloneloadServeConstraint[bendersIter],
			cloneloadServeReactiveConstraint[bendersIter],
			cloneswitchImplicationConstraint1[bendersIter],
			//cloneswitchImplicationConstraint2[bendersIter],
			clonehardenConstraint[bendersIter],
			clonevoltageConstraint[bendersIter],
			clonevoltageOffsetConstraint[bendersIter],
			clonevoltageVariationConstraint1[bendersIter],
			clonevoltageVariationConstraint2[bendersIter],
			clonelineUseVariable[bendersIter], // VARIABLES
			clonelineExistsVariable[bendersIter], // VARIABLES
			clonelineDirectionVariable[bendersIter],
			clonelineCycleVariable[bendersIter],
			clonelineHardenVariable[bendersIter],
			cloneswitchUseVariable[bendersIter],
			cloneswitchCycleVariable[bendersIter],
			cloneflowVariable[bendersIter],
			cloneflowReactiveVariable[bendersIter],
			clonevoltageVariable[bendersIter],
			clonevoltageOffsetVariable[bendersIter],
			cloneloadVariable[bendersIter],
			cloneloadReactiveVariable[bendersIter],
			cloneloadServeVariable[bendersIter],
			clonegeneratorVariable[bendersIter],
			clonegeneratorReactiveVariable[bendersIter],
			clonefacilityVariable[bendersIter],
			numDamagedEdges,
			NumPhases,
			ScenarioIndex,
			NumScenarios,
			discreteMicrogrid,
			damage_str,
			root_str,
			HardenedDisabled,
			DISABLED[ScenarioIndex],
			HARDENED_DISABLED[ScenarioIndex],
			newImpedanceMultiplier,
			LDFIndicator);


	
		SubProblem[NumScenarios].add(clonelineExistsConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonelineConstraint1[bendersIter]);
		SubProblem[NumScenarios].add(clonelineConstraint2[bendersIter]);
		SubProblem[NumScenarios].add(clonelineReactiveConstraint1[bendersIter]);
		SubProblem[NumScenarios].add(clonelineReactiveConstraint2[bendersIter]);
		SubProblem[NumScenarios].add(clonedirectionConstraint[bendersIter]);
		SubProblem[NumScenarios].add(cloneswitchConstraint1[bendersIter]);
		SubProblem[NumScenarios].add(cloneswitchConstraint2[bendersIter]);
		SubProblem[NumScenarios].add(cloneswitchReactiveConstraint1[bendersIter]);
		SubProblem[NumScenarios].add(cloneswitchReactiveConstraint2[bendersIter]);
		SubProblem[NumScenarios].add(clonevariationConstraint1[bendersIter]);
		SubProblem[NumScenarios].add(clonevariationConstraint2[bendersIter]);
		SubProblem[NumScenarios].add(clonevariationReactiveConstraint1[bendersIter]);
		SubProblem[NumScenarios].add(clonevariationReactiveConstraint2[bendersIter]);
		SubProblem[NumScenarios].add(clonedamageConstraint[bendersIter]);
		SubProblem[NumScenarios].add(cloneloadConstraint[bendersIter]);
		SubProblem[NumScenarios].add(cloneloadReactiveConstraint[bendersIter]);
		//SubProblem[NumScenarios].add(clonegenerationConstraint[bendersIter]);
		//SubProblem[NumScenarios].add(clonegenerationReactiveConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonebalanceConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonebalanceReactiveConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonemicrogridConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonecycleConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonecycleSubConstraint[bendersIter]);
		SubProblem[NumScenarios].add(cloneswitchableConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonecriticalServeConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonecriticalServeReactiveConstraint[bendersIter]);
		SubProblem[NumScenarios].add(cloneloadServeConstraint[bendersIter]);
		SubProblem[NumScenarios].add(cloneloadServeReactiveConstraint[bendersIter]);
		SubProblem[NumScenarios].add(cloneswitchImplicationConstraint1[bendersIter]);
		//SubProblem[NumScenarios].add(cloneswitchImplicationConstraint2[bendersIter]);
		SubProblem[NumScenarios].add(clonehardenConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonevoltageConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonevoltageOffsetConstraint[bendersIter]);
		SubProblem[NumScenarios].add(clonevoltageVariationConstraint1[bendersIter]);
		SubProblem[NumScenarios].add(clonevoltageVariationConstraint2[bendersIter]);

		SubProblem[NumScenarios].add(clonelineUseVariable[bendersIter]);
		SubProblem[NumScenarios].add(clonelineExistsVariable[bendersIter]);
		SubProblem[NumScenarios].add(clonelineDirectionVariable[bendersIter]);
		SubProblem[NumScenarios].add(clonelineCycleVariable[bendersIter]);
		SubProblem[NumScenarios].add(clonelineHardenVariable[bendersIter]);
		SubProblem[NumScenarios].add(cloneswitchUseVariable[bendersIter]);
		SubProblem[NumScenarios].add(cloneswitchCycleVariable[bendersIter]);
		SubProblem[NumScenarios].add(cloneflowVariable[bendersIter]);
		SubProblem[NumScenarios].add(cloneflowReactiveVariable[bendersIter]);
		SubProblem[NumScenarios].add(clonevoltageVariable[bendersIter]);
		SubProblem[NumScenarios].add(clonevoltageOffsetVariable[bendersIter]);
		SubProblem[NumScenarios].add(cloneloadVariable[bendersIter]);
		SubProblem[NumScenarios].add(cloneloadReactiveVariable[bendersIter]);
		SubProblem[NumScenarios].add(cloneloadServeVariable[bendersIter]);
		SubProblem[NumScenarios].add(clonegeneratorVariable[bendersIter]);
		SubProblem[NumScenarios].add(clonegeneratorReactiveVariable[bendersIter]);
		//SubProblem[NumScenarios].add(clonemicrogridVariable[bendersIter]);
		//SubProblem[NumScenarios].add(clonefacilityVariable[bendersIter]);



		//cout << "LINES!" << endl;
		for (int j = 0; j < NumEdges; ++j) {
			BaseLineVarConstraint.add(IloRange(env, 0.0, IloInfinity));
			BaseLineVarConstraint[bendersIter*NumEdges + j].setLinearCoef(lineUseVariable[NumScenarios][j], 1.0);
			BaseLineVarConstraint[bendersIter*NumEdges + j].setLinearCoef(clonelineUseVariable[bendersIter][j], -1.0);
			SubProblem[NumScenarios].add(BaseLineVarConstraint[bendersIter*NumEdges + j]);
		}
		//cout << "SWITCHES!" << endl;
		for (int j = 0; j < NumEdges; ++j) {
			BaseSwitchVarConstraint.add(IloRange(env, 0.0, IloInfinity));
			BaseSwitchVarConstraint[bendersIter*NumEdges + j].setLinearCoef(switchUseVariable[NumScenarios][j], 1.0);
			BaseSwitchVarConstraint[bendersIter*NumEdges + j].setLinearCoef(cloneswitchUseVariable[bendersIter][j], -1.0);
			SubProblem[NumScenarios].add(BaseSwitchVarConstraint[bendersIter*NumEdges + j]);
		}
		//cout << "HARDENS!" << endl;
		for (int j = 0; j < NumEdges; ++j) {
			BaseHardenVarConstraint.add(IloRange(env, 0.0, IloInfinity));
			BaseHardenVarConstraint[bendersIter*NumEdges + j].setLinearCoef(lineHardenVariable[NumScenarios][j], 1.0);
			BaseHardenVarConstraint[bendersIter*NumEdges + j].setLinearCoef(clonelineHardenVariable[bendersIter][j], -1.0);
			SubProblem[NumScenarios].add(BaseHardenVarConstraint[bendersIter*NumEdges + j]);
		}
		//cout << "FACILITIES!" << endl;
		int numCon = -1;
		int numPrevCon = BaseFacilityVarConstraint.getSize();
		for (int j = 0; j < NODES.size(); ++j) {
			for (int k = 0; k < NumPhases; ++k) {
				++numCon;
				BaseFacilityVarConstraint.add(IloRange(env, -MAXREALPHASE[j][k], IloInfinity));
				BaseFacilityVarConstraint[numPrevCon + numCon].setLinearCoef(clonegeneratorVariable[bendersIter][j*NumPhases + k], -1.0);
				for (int l = 0; l < GENERATORS.size(); ++l) {	
					if (hashTableVertex[GENERATORS[l].node_id] == j && GENERATORS[l].hasphase[k]) {
						BaseFacilityVarConstraint[numPrevCon + numCon].setLinearCoef(facilityVariable[NumScenarios][l], MAX_ADDED[l]);
					}
				}
				SubProblem[NumScenarios].add(BaseFacilityVarConstraint[numPrevCon + numCon]);	
			}
		}
	}



	// -- GET LP RELAXATION ----------------------------------------------------------------
	vector<double> LPSolution(ColSize, 0.0);

	solveRelaxed(
		SubProblem[NumScenarios],
		SubSolver[NumScenarios],
		lineUseVariable[NumScenarios], // VARIABLES
		lineExistsVariable[NumScenarios], // VARIABLES
		lineDirectionVariable[NumScenarios],
		lineCycleVariable[NumScenarios],
		lineHardenVariable[NumScenarios],
		switchUseVariable[NumScenarios],
		switchCycleVariable[NumScenarios],
		flowVariable[NumScenarios],
		flowReactiveVariable[NumScenarios],
		voltageVariable[NumScenarios],
		voltageOffsetVariable[NumScenarios],
		loadVariable[NumScenarios],
		loadReactiveVariable[NumScenarios],
		loadServeVariable[NumScenarios],
		generatorVariable[NumScenarios],
		generatorReactiveVariable[NumScenarios],
		facilityVariable[NumScenarios],
		clonelineUseVariable, // SUB PROBLEM VARIABLES
		clonelineExistsVariable, // SUB PROBLEM VARIABLES
		clonelineDirectionVariable,
		clonelineCycleVariable,
		clonelineHardenVariable,
		cloneswitchUseVariable,
		cloneswitchCycleVariable,
		cloneflowVariable,
		cloneflowReactiveVariable,
		clonevoltageVariable,
		clonevoltageOffsetVariable,
		cloneloadVariable,
		cloneloadReactiveVariable,
		cloneloadServeVariable,
		clonegeneratorVariable,
		clonegeneratorReactiveVariable,
		LPSolution,
		NumScenarios,
		NumEdges,
		NumGenerators,
		NumPhases);
	// -------------------------------------------------------------------------------------


	vector<double> incumbentSolution(ColSize, 0.0);
	for (int j = 0; j < ColSize; ++j) {
		incumbentSolution[j] = MIN_FEAS_SOL[j];
	}

	vector<IloNumArray> cloneLineExistsVals(0);
	vector<IloNumArray> cloneLineUseVals(0);
	vector<IloNumArray> cloneLineDirectionVals(0);
	vector<IloNumArray> cloneSwitchUseVals(0);
	vector<IloNumArray> cloneGeneratorVals(0);
	vector<IloNumArray> cloneGeneratorReactiveVals(0);
	vector<IloNumArray> cloneLoadVals(0);
	vector<IloNumArray> cloneLoadReactiveVals(0);
	vector<IloNumArray> cloneLoadServeVals(0);
	vector<IloNumArray> cloneVoltageVals(0);

	time_t begin = clock();
	VariableNeighborhoodSearchSolver(
		MIN_FEAS_OBJ,
		incumbentSolution,
		LPSolution,
		SubProblem,
		SubSolver,
		hashTableVertex,
		hashTableEdge,
		hashTableGenerator,
		hashTableLoad,
		NODES,
		EDGES,
		GENERATORS,
		LOADS,
		LINECODES,
		LINE_CONSTRUCTION_COST,
		MICROGRID_COST,
		MICROGRID_FIXED_COST,
		IS_CRITICAL_LOAD,
		MAX_MICROGRID,
		HARDEN_COST,
		LINE_SWITCH_COST,
		LoadMet,
		CriticalLoadMet,
		PhaseVariation,
		G,
		CYCLES,
		lineUseVariable, // VARIABLES
		lineExistsVariable, // VARIABLES
		lineDirectionVariable,
		lineCycleVariable,
		lineHardenVariable,
		switchUseVariable,
		switchCycleVariable,
		flowVariable,
		flowReactiveVariable,
		voltageVariable,
		voltageOffsetVariable,
		loadVariable,
		loadReactiveVariable,
		loadServeVariable,
		generatorVariable,
		generatorReactiveVariable,
		facilityVariable,
		clonelineUseVariable,
		clonelineExistsVariable,
		clonelineDirectionVariable,
		cloneswitchUseVariable,
		clonegeneratorVariable,
		clonegeneratorReactiveVariable,
		clonevoltageVariable,
		clonevoltageOffsetVariable,
		cloneloadVariable,
		cloneloadReactiveVariable,
		cloneloadServeVariable,
		cloneLineExistsVals,
		cloneLineUseVals,
		cloneLineDirectionVals,
		cloneSwitchUseVals,
		cloneGeneratorVals,
		cloneGeneratorReactiveVals,
		cloneLoadVals,
		cloneLoadReactiveVals,
		cloneLoadServeVals,
		cloneVoltageVals,
		NumScenarios,
		NumDamagedEdges,
		NumPhases,
		NumEdges,
		NumGenerators,
		NumScenarios,
		discreteMicrogrid,
		damage_str,
		extent_str,
		root_str,
		HardenedDisabled,
		48.0 * 60.0 * 60.0);
	time_t end = clock();

	double incumbentObjective = std::inner_product(ObjectiveVector.begin(), ObjectiveVector.end(), incumbentSolution.begin(), 0.0);
	if (CHECK_RESULTS) {
		env = SubSolver[NumScenarios].cplex.getEnv();
		IloNumArray lineusevals(env);
		IloNumArray lineswitchvals(env);
		IloNumArray linehardenvals(env);
		//IloNumArray microgridvals(env);				
		IloNumArray facilityvals(env);

		for (int j = 0; j < NumEdges; ++j) {
			lineusevals.add(incumbentSolution[j]);
			lineswitchvals.add(incumbentSolution[j + NumEdges]);
			linehardenvals.add(incumbentSolution[j + 2*NumEdges]);
		}
		for (int j = 0; j < NumGenerators; ++j) {
			/*for (int k = 0; k < NumPhases; ++k) {
				microgridvals.add(incumbentSolution[j*NumPhases + k + 3*NumEdges]);
			}*/
			facilityvals.add(incumbentSolution[j + 3*NumEdges]);
		}

		std::stringstream ss1;
		ss1 << CriticalLoadMet;
		std::stringstream ss2;
		ss2 << discreteMicrogrid;
		string critical_str = ss1.str() + "_" + ss2.str();

		string str;
		str = "results/VNS/" + root_str + "solution_" + critical_str + "_" + damage_str + ".txt";
		
		cout << str << endl;
		ofstream outputSolution(str.c_str());

		outputSolution << "FLOW VALUES: " << endl;
		for (int j = 0; j < EDGES.size(); ++j) {
			outputSolution << "EDGE " << EDGES[j].id << " (" << lineusevals[j] << "," << lineswitchvals[j] << "," << linehardenvals[j] << ")" << endl;
		}
		outputSolution << "GENERATION VALUES: " << endl;
		for (int j = 0; j < GENERATORS.size(); ++j) {
			outputSolution << "GENERATOR " << GENERATORS[j].id << " (" << facilityvals[j] << ")" << endl;
			
		}

		/*for (int j = 0; j < EDGES.size(); ++j) {
			if (lineusevals[j] == 1) {
				cout << EDGES[j].node1id << " " << EDGES[j].node2id << endl;
			}
		}*/

		outputSolution << "STATUS: Feasible OBJECTIVE: " << incumbentObjective << endl;
		outputSolution << "CPUTIME: " << (double) (end - begin) / (double) CLOCKS_PER_SEC << " NUM SCENARIOS GENERATED: " << bendersIter+1 << endl;


		// CHECK IF THERE ARE CYCLES IN THE SOLUTION:
		bool CHECK_CYCLES = true;
		if (CHECK_CYCLES) {
			vector<bool> edgeIncluded (EDGES.size(), false);
			for (int j = 0; j < EDGES.size(); ++j) {
				if (fabs(lineusevals[j] - 1.0) < 1e-4 && fabs(lineswitchvals[j] - 0.0) < 1e-4) {
					edgeIncluded[j] = true;
				}
			}
			// CREATE GRAPH
			vector<graph_vertex> GSub(NODES.size());
			for (int j = 0; j < EDGES.size(); ++j) {
				if (edgeIncluded[j]) {
					GSub[hashTableVertex[EDGES[j].node1id]].AdjList.push_back(hashTableVertex[EDGES[j].node2id]);
					GSub[hashTableVertex[EDGES[j].node1id]].EdgeID.push_back(EDGES[j].id);
					GSub[hashTableVertex[EDGES[j].node2id]].AdjList.push_back(hashTableVertex[EDGES[j].node1id]);
					GSub[hashTableVertex[EDGES[j].node2id]].EdgeID.push_back(EDGES[j].id);
				}
			}

			// CREATE CYCLE LIST AND FIND ALL CYCLES
			vector<vector<int> > CYCLESSub;
			detectCycles(GSub, CYCLESSub);
			outputSolution << "SUBGRAPH HAS " << CYCLESSub.size() << " CYCLES" << endl;
			outputSolution << "MINIMUM FEASIBLE OBJECTIVE: " << MIN_FEAS_OBJ << endl;
		}

	}


}

void readGridLabDOutput(
	vector<nodeData> &NODES,
	vector<edgeData> &EDGES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	map<int, lineCodeData> &LINECODES,
	map<string, int> &hashTableVertex,
	map<string, int> &hashTableEdge,
	map<string, int> &hashTableGenerator,
	map<string, int> &hashTableLoad,
	vector<int> &scenarioList,
	vector<list<string> > &nodesList,
	vector<string> &infNodeList,
	vector<list<string> > &edgesList,
	int NumScenarios,
	int NumNodes,
	int NumEdges,
	int NumPhases,
	vector<vector<int> > &allForests,
	vector<vector<bool> > &activeGenerators,
	string gridlabd_path) {

	// CUT LIST:
	int numPrevComp = scenarioList.size();

	//cout << "STARTED READING GLD OUTPUT" << endl;
        string gld_output_path = gridlabd_path + string("/gld_output");
	ifstream file(gld_output_path.c_str());
	string line;
	while (!file.eof()) {
		getline(file, line);
		istringstream ss(line);
		string s;
		getline( ss, s, '=');
		getline( ss, s, '\n');
		if (atoi(s.c_str()) < 1) break;
		int ind = atoi(s.c_str()) - 1;
		//cout << "Scenario #" << ind << endl;
		getline(file, line);
		istringstream ssNodes(line);
		list<string> nodes; 
		while (getline(ssNodes, s, ' ')) {
			if (s == "10000") 
				nodes.push_back("sourcebus");
			else
				nodes.push_back(s);
			//cout << s << endl;
		}
		getline(file, line);
		istringstream ssinfNodes(line);
		list<string> infNodes;
		while (getline(ssinfNodes, s, ' ')) {
			if (s == "10000") 
				infNodes.push_back("sourcebus");
			else
				infNodes.push_back(s);
			//cout << s << endl;
		}
		//cout << "INFEASIBLE " << infNode << endl;	

		for (list<string>::iterator itMain = infNodes.begin();
			itMain != infNodes.end(); ++itMain) {
			list<string>::iterator itBegin = find(nodes.begin(), nodes.end(), *itMain);

			list<string> nodesSub;
			for (list<string>::iterator itSub = nodes.begin();
				itSub != itBegin; ++itSub) {
                                if (find(infNodes.begin(), infNodes.end(), *itSub) == infNodes.end()) {
					nodesSub.push_back(*itSub);	
                                }
			}
			nodesSub.push_back(*itBegin);

			cout << "NEW INFEASIBLE SUBSET: ";
			for (list<string>::iterator myit = nodesSub.begin(); myit != nodesSub.end(); ++myit) {
				cout << *myit << " ";
			}
			cout << endl;

			// ADD GENERATORS TO THE SUBSETS:
			for (int l = 0; l < NumNodes; ++l) {
				//cout << activeGenerators[ind][l] << " ";
				if (activeGenerators[ind][l]) {// && find(nodes.begin(), nodes.end(), NODES[l].id) == nodes.end()) {
					nodesSub.remove(NODES[l].id);
					nodesSub.push_front(NODES[l].id);
				}
			}
			//cout << endl;

                        //if (find(nodesList.begin(), nodesList.end(), nodesSub) == nodesList.end()) {
				scenarioList.push_back(ind);
				nodesList.push_back(nodesSub);
				infNodeList.push_back(*itBegin);
			//}
                }


	}
	file.close();	
	// COMPUTE SUB COMPONENTS FOR EDGES

	//cout << "GLD OUTPUT DONE" << endl;	
	for (int l = numPrevComp; l < scenarioList.size(); ++l) {
		// CREATE GRAPH
		int i = scenarioList[l];
		vector<graph_vertex> GSub(NumNodes);
		for (int j = 0; j < EDGES.size(); ++j) {
			if (allForests[i][j] != 0) {
				GSub[hashTableVertex[EDGES[j].node1id]].AdjList.push_back(hashTableVertex[EDGES[j].node2id]);
				GSub[hashTableVertex[EDGES[j].node1id]].EdgeID.push_back(EDGES[j].id);
				GSub[hashTableVertex[EDGES[j].node2id]].AdjList.push_back(hashTableVertex[EDGES[j].node1id]);
				GSub[hashTableVertex[EDGES[j].node2id]].EdgeID.push_back(EDGES[j].id);
			}
		}
		int source = hashTableVertex[*(nodesList[l].begin())];
		list<vector<int> > currIndices;
		list<string> currList;
	
		//cout << "MINIMAL INFEASIBLE LIST" << endl;	
		// CREATE MINIMAL INFEASIBLE SET OF NODES
		for (list<string>::iterator it = nodesList[l].begin();
			it != nodesList[l].end(); ++it) {
			//cout << *it << " ";
			GSub[hashTableVertex[*it]].included = true;
			/*if (*it == infNodeList[l]) {
				break;
			}*/
		}
		//cout << endl;

		//cout << "COMPUTE SUBCOMPONENT" << endl;
		findSubComponent(GSub, source, currIndices);
		for (int j = 0; j < EDGES.size(); ++j) {
			vector<int> edge(2);
			edge[0] = hashTableVertex[EDGES[j].node1id];
			edge[1] = hashTableVertex[EDGES[j].node2id];
			std::sort(edge.begin(), edge.end());
			if (find(currIndices.begin(), currIndices.end(), edge) != currIndices.end() && allForests[scenarioList[l]][j] != 0) {
				currList.push_back(EDGES[j].id);
				//cout << EDGES[j].id << " ";
			}
		}
		//cout << endl;
		edgesList.push_back(currList);
		//cout << "COMPONENT DONE" << endl;
	}



}


void populateGridLabDInput(
	vector<nodeData> &NODES,
	vector<edgeData> &EDGES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	map<int, lineCodeData> &LINECODES,
	map<string, int> &hashTableVertex,
	map<string, int> &hashTableEdge,
	map<string, int> &hashTableGenerator,
	map<string, int> &hashTableLoad,
	vector<IloNumArray> &cloneLineUseVals,
	vector<IloNumArray> &cloneLineDirectionVals,
	vector<IloNumArray> &cloneSwitchUseVals,
	vector<IloNumArray> &cloneGeneratorVals,
	vector<IloNumArray> &cloneGeneratorReactiveVals,
	vector<IloNumArray> &cloneLoadVals,
	vector<IloNumArray> &cloneLoadReactiveVals,
	vector<IloNumArray> &cloneLoadServeVals,
	vector<IloNumArray> &cloneVoltageVals,
	vector<IloNumArray> &subLineUseVals,
	vector<IloNumArray> &subLineDirectionVals,
	vector<IloNumArray> &subSwitchUseVals,
	vector<IloNumArray> &subGeneratorVals,
	vector<IloNumArray> &subGeneratorReactiveVals,
	vector<IloNumArray> &subLoadVals,
	vector<IloNumArray> &subLoadReactiveVals,
	vector<IloNumArray> &subLoadServeVals,
	vector<IloNumArray> &subVoltageVals,
	vector<bool> &ScenarioIncluded,
	vector<int> &addedScenarios,
	int NumScenarios,
	int NumNodes,
	int NumEdges,
	int NumLoads,
	int NumPhases,
	vector<vector<int> > &allForests,
	vector<vector<int> > &allLoads,
	vector<vector<vector<double> > > &allVoltage,
	vector<vector<vector<double> > > &allRealGen,
	vector<vector<vector<double> > > &allReactiveGen,
	vector<vector<vector<double> > > &allRealLoad,
	vector<vector<vector<double> > > &allReactiveLoad,
	vector<vector<bool> > &activeGenerators,
        vector<vector<bool> > &activeLoads) {
	// FIND ALL FORESTS:
	for (int j = 0; j < cloneLineUseVals.size(); ++j) {
		for (int e = 0; e < NumEdges; ++e) {
			if (fabs((double)cloneLineUseVals[j][e] - 1.0) < 1e-4 && fabs((double)cloneSwitchUseVals[j][e] - 0.0) < 1e-4) {
				if (fabs((double) cloneLineDirectionVals[j][2*e] - 1.0) < 1e-4) {
					allForests[addedScenarios[j]][e] = -1;
				} else if (fabs((double) cloneLineDirectionVals[j][2*e+1] - 1.0) < 1e-4) {
					allForests[addedScenarios[j]][e] = 1;

				}
				//cout << "1 ";
			} else {
				//cout << "0 ";
			}
		}
		//cout << endl;
		for (int l = 0; l < NumNodes; ++l) {
			for (int k = 0; k < NumPhases; ++k) {
				//cout << cloneVoltageVals[j].getSize() << " " << (l*NumPhases+k) << endl;
				//cout << addedScenarios[j] << endl;
				//cout << pow(cloneVoltageVals[j][l*NumPhases + k]/1e+3,0.5) << endl;
				//cout << NODES[l].refVoltage << endl;
				allVoltage[addedScenarios[j]][l][k] = pow(cloneVoltageVals[j][l*NumPhases + k]/1e+3,0.5) / NODES[l].refVoltage;
				//cout << "AHOY!!" << endl;
				allRealGen[addedScenarios[j]][l][k] = cloneGeneratorVals[j][l*NumPhases + k];
				//cout << "AHOY!!" << endl;
				allReactiveGen[addedScenarios[j]][l][k] = cloneGeneratorReactiveVals[j][l*NumPhases + k];
				//cout << "AHOY!!" << endl;
				allRealLoad[addedScenarios[j]][l][k] = cloneLoadVals[j][l*NumPhases + k];
				//cout << "AHOY!!" << endl;
				allReactiveLoad[addedScenarios[j]][l][k] = cloneLoadReactiveVals[j][l*NumPhases + k];
				//cout << "AHOY!!" << endl;
			}
		}
 		for (int e = 0; e < NumLoads; ++e) {
			if (fabs((double)cloneLoadServeVals[j][e] - 1.0) < 1e-4) {
				allLoads[addedScenarios[j]][e] = 1;
			} else {
				//cout << "0 ";
			}
		}
	}
	for (int j = 0; j < NumScenarios; ++j) {
		if (!ScenarioIncluded[j]) {
			for (int e = 0; e < NumEdges; ++e) {
				if (fabs((double)subLineUseVals[j][e] - 1.0) < 1e-4 && fabs((double)subSwitchUseVals[j][e] - 0.0) < 1e-4) {
					if (fabs((double) subLineDirectionVals[j][2*e] - 1.0) < 1e-4) {
						allForests[j][e] = -1;
					} else if (fabs((double) subLineDirectionVals[j][2*e+1] - 1.0) < 1e-4) {
						allForests[j][e] = 1;
					}
				}
			} 
			for (int l = 0; l < NumNodes; ++l) {
				for (int k = 0; k < NumPhases; ++k) {
					allVoltage[j][l][k] = pow(subVoltageVals[j][l*NumPhases + k]/1e+3,0.5) / NODES[l].refVoltage;
					allRealGen[j][l][k] = subGeneratorVals[j][l*NumPhases + k];
					allReactiveGen[j][l][k] = subGeneratorReactiveVals[j][l*NumPhases + k];
					allRealLoad[j][l][k] = subLoadVals[j][l*NumPhases + k];
					allReactiveLoad[j][l][k] = subLoadReactiveVals[j][l*NumPhases + k];
				}
			}
			for (int e = 0; e < NumLoads; ++e) {
				if (fabs((double)subLoadServeVals[j][e] - 1.0) < 1e-4) {
					allLoads[j][e] = 1;
				}
			}
		}
	}
	for (int j = 0; j < NumScenarios; ++j) {
		for (int l = 0; l < NumNodes; ++l) {
			for (int k = 0; k < NumPhases; ++k) {
				if (allRealGen[j][l][k] > 1e-4 || allReactiveGen[j][l][k] > 1e-4) {
					activeGenerators[j][l] = true;
				}
				if (allRealLoad[j][l][k] > 1e-4 || allReactiveLoad[j][l][k] > 1e-4) {
					activeLoads[j][l] = true;
				}

			}
		}
	}


}

void GridLabDInput(
	vector<nodeData> &NODES,
	vector<edgeData> &EDGES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	map<int, lineCodeData> &LINECODES,
	map<string, int> &hashTableVertex,
	map<string, int> &hashTableEdge,
	map<string, int> &hashTableGenerator,
	map<string, int> &hashTableLoad,
	vector<bool> &ScenarioIncluded,
	vector<int> &addedScenarios,
	int NumScenarios,
	int NumNodes,
	int NumEdges,
	int NumPhases,
	vector<vector<int> > &allForests,
	vector<vector<vector<double> > > &allVoltage,
	vector<vector<vector<double> > > &allRealGen,
	vector<vector<vector<double> > > &allReactiveGen,
	vector<vector<vector<double> > > &allRealLoad,
	vector<vector<vector<double> > > &allReactiveLoad,
	int ComponentOrder,
	int gridlabdindex,
	string gridlabd_path) {

	cout << "GRIDLABD OUTPUT" << endl;


	
	cout << "NUMBER OF DISTINCT TREES: " << allForests.size() << endl;
        stringstream ss;
        ss << gridlabdindex;
        std::string output_file = gridlabd_path + string("/forests_") + ss.str() + string(".txt");
	ofstream out1(output_file.c_str());
	for (int e = 0; e < NumEdges; ++e) {
		out1 << EDGES[e].id << ",";
		for (int i = 0; i < allForests.size(); ++i) {
			int edgeVal = allForests[i][e]; 
			int suppEdgeVal = true;
			string node1id = EDGES[e].node1id;
			string node2id = EDGES[e].node2id;
			if (node1id[node1id.size()-1] == 'r') { 
				for (int j = 0; j < NumEdges; ++j) {
					if (e != j && (EDGES[j].node1id == node1id ||
							EDGES[j].node2id == node1id)) {
						if (allForests[i][j] != 0) {
							suppEdgeVal = false;
						}
					}
				}
				if (suppEdgeVal) {
					edgeVal = 0;
				}
			} else if (node2id[node2id.size()-1] == 'r') { 
				for (int j = 0; j < NumEdges; ++j) {
					if (e != j && (EDGES[j].node1id == node2id ||
							EDGES[j].node2id == node2id)) {
						if (allForests[i][j] != 0) {
							suppEdgeVal = false;
						}
					}
				}
				if (suppEdgeVal) {
					edgeVal = 0;
				}
			}
			out1 << edgeVal << ",";
		}
		out1 << endl;
	}
	out1.close();


	vector<int> FirstComponent;
	vector<int> FirstEdges;


	/*ofstream out("COMPONENTS.txt", ios::app);

	out << "NEW ITERATION" << endl;

	// CHECK IF THERE ARE CYCLES IN THE SOLUTION:
	for (int i = 0; i < NumScenarios; ++i) {
		// CREATE GRAPH
		vector<graph_vertex> GSub(NODES.size());
		for (int j = 0; j < EDGES.size(); ++j) {
			if (allForests[i][j] != 0) {
				GSub[hashTableVertex[EDGES[j].node1id]].AdjList.push_back(hashTableVertex[EDGES[j].node2id]);
				GSub[hashTableVertex[EDGES[j].node1id]].EdgeID.push_back(EDGES[j].id);
				GSub[hashTableVertex[EDGES[j].node2id]].AdjList.push_back(hashTableVertex[EDGES[j].node1id]);
				GSub[hashTableVertex[EDGES[j].node2id]].EdgeID.push_back(EDGES[j].id);
			}
		}

		// CREATE CYCLE LIST AND FIND ALL CYCLES
		vector<vector<int> > CYCLESSub;
		detectCycles(GSub, CYCLESSub);
		cout << "SUBGRAPH " << i << " HAS " << CYCLESSub.size() << " CYCLES" << endl;
		list<vector<int> > components;
		findComponents(GSub, components);

		for (list<vector<int> >::iterator it = components.begin();
			it != components.end(); ++it) {
			int numGenerators = 0;

			vector<double> componentRealGen(NumPhases, 0.0);				
			vector<double> componentReactiveLoad(NumPhases, 0.0);
			if (it == components.begin()) cout << "FIRST COMPONENT: SIZE = " << it->size() << ", ";
			for (int j = 0; j < it->size(); ++j) {
				for (int k = 0; k < NumPhases; ++k) {
					componentRealGen[k] += allRealGen[i][(*it)[j]][k];
					componentReactiveLoad[k] += allReactiveLoad[i][(*it)[j]][k];
				}		
			}
			bool activeComponent = false;
			for (int k = 0; k < NumPhases; ++k) {
				if (componentRealGen[k] > 1e-4) {
					activeComponent = true;
				}
			}
			if (activeComponent) {
				out << "NEW ACTIVE COMPONENT: " << endl << "NODES: ";
				for (int j = 0; j < it->size(); ++j) {
					FirstComponent.push_back((*it)[j]);
					out << NODES[(*it)[j]].id << " ";
				}
				out << endl << "EDGES: ";
				for (int j = 0; j < NumEdges; ++j) {
					int node1ind = hashTableVertex[EDGES[j].node1id];	
					int node2ind = hashTableVertex[EDGES[j].node2id];
					bool cond1 = std::find(it->begin(), it->end(), node1ind) != it->end();	
					bool cond2 = std::find(it->begin(), it->end(), node2ind) != it->end();	
					if (cond1 && cond2) {
						//if (allForests[0][j] != 0) {
							FirstEdges.push_back(j);
							out << EDGES[j].id << " ";
						//}
					}
				}
				out << endl;
			}*/
			/*for (int j = 0; j < it->size(); ++j) {
				for (int k = 0; k < NumPhases; ++k) {
					if (componentRealGen[k] > 1e-4) {
						allReactiveGen[i][(*it)[j]][k] = (allRealGen[i][(*it)[j]][k] / componentRealGen[k]) * componentReactiveLoad[k];
					}
				}
			}*/
		/*}
	}
	out.close();*/

	/*std::sort(FirstComponent.begin(), FirstComponent.end());
	std::sort(FirstEdges.begin(), FirstEdges.end());

	map<vector<int>, set<vector<int> > >::iterator it = componentSolutions.find(FirstComponent);

	if (it == componentSolutions.end()) {
		set<vector<int> > myset;
		myset.insert(FirstEdges);
		componentSolutions.insert(pair<vector<int>, set<vector<int> > >(FirstComponent, myset));
	} else {
		if (it->second.find(FirstEdges) != it->second.end()) {
			cout << "********* IMPOSSIBLE PREVIOUS COMPONENT **********" << endl;
		} else {
			it->second.insert(FirstEdges);
		}
	}*/

	int NumComponents = 3;
	vector<sortWithIndices> voltageSumPre(NumNodes);
	vector<sortWithIndices> voltageCompSum(NumComponents);
	vector<sortWithIndices> voltageSum(NumNodes);
	

	if (ComponentOrder == 0) {
		for (int j = 0; j < NumNodes; ++j) {
			voltageSumPre[j].ind = j;
			voltageSumPre[j].val = 0.0;
			voltageSumPre[j].tie = 0.0;
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					voltageSumPre[j].val += allVoltage[i][j][k];
				}
			}
		}
		for (int j = 0; j < NumComponents; ++j) {	
			std::sort(voltageSumPre.begin()+(1+NumNodes/3*j), voltageSumPre.begin()+(1+NumNodes/3*(j+1)), mycomperator);
			voltageCompSum[j].ind = j;
			voltageCompSum[j].val = 0.0;
			voltageCompSum[j].tie = 0.0;
			for (int i = 0; i < NumScenarios; ++i) {
				for (int l = 0; l < (NumNodes-1)/3; ++l) {
					for (int k = 0; k < NumPhases; ++k) {
						voltageCompSum[j].val += allVoltage[i][1+NumNodes/3*j+l][k];
					}
				}
			}
		}
		std::sort(voltageCompSum.begin(), voltageCompSum.end(), mycomperator);
		voltageSum[0].ind = voltageSumPre[0].ind;
		voltageSum[0].val = voltageSumPre[0].val;
		voltageSum[0].tie = voltageSumPre[0].tie;
		for (int j = 0; j < NumComponents; ++j) {
			for (int l = 0; l < (NumNodes-1)/3; ++l) {
				voltageSum[1+NumNodes/3*j+l].ind = voltageSumPre[1+NumNodes/3*(voltageCompSum[j].ind)+l].ind;
				voltageSum[1+NumNodes/3*j+l].val = voltageSumPre[1+NumNodes/3*(voltageCompSum[j].ind)+l].val;
				voltageSum[1+NumNodes/3*j+l].tie = voltageSumPre[1+NumNodes/3*(voltageCompSum[j].ind)+l].tie;
			}
		}
	} else if (ComponentOrder == 1) {
		for (int j = 0; j < NumNodes; ++j) {
			voltageSum[j].ind = j;
			voltageSum[j].val = 0.0;
			voltageSum[j].tie = 0.0;
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					voltageSum[j].val += allVoltage[i][j][k];
				}
			}
		}
		std::sort(voltageSum.begin()+1, voltageSum.end(), mycomperator);
	} else if (ComponentOrder == 2) {
		for (int j = 0; j < NumNodes; ++j) {
			voltageSum[j].ind = j;
			voltageSum[j].val = 0.0;
			voltageSum[j].tie = 0.0;
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					voltageSum[j].val += allVoltage[i][j][k];
				}
			}
		}
		std::random_shuffle(voltageSum.begin()+1, voltageSum.end());
	}

	


	output_file = gridlabd_path + string("/voltages_") + ss.str() + string(".txt");
	ofstream out7(output_file.c_str());
	for (int j = 0; j < NumNodes; ++j) {
		int e = voltageSum[j].ind;
		if (find(NODES[e].id.begin(), NODES[e].id.end(), 'r') == NODES[e].id.end() || NODES[e].id[0] == 's') {
			out7 << NODES[e].id << ",";
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					out7 << allVoltage[i][e][k] << ",";
				}
			}
			out7 << endl;
		}
	}
	out7.close();
	output_file = gridlabd_path + string("/loadsReal_") + ss.str() + string(".txt");
	ofstream out2(output_file.c_str());
	for (int j = 0; j < NumNodes; ++j) {
		int e = voltageSum[j].ind;
		if (find(NODES[e].id.begin(), NODES[e].id.end(), 'r') == NODES[e].id.end() || NODES[e].id[0] == 's') {
			out2 << NODES[e].id << ",";
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					out2 << allRealLoad[i][e][k] << ",";
				}
			}
			out2 << endl;
		}
	}
	out2.close();
	output_file = gridlabd_path + string("/loadsReactive_") + ss.str() + string(".txt");
	ofstream out3(output_file.c_str());
	for (int j = 0; j < NumNodes; ++j) {
		int e = voltageSum[j].ind;
		if (find(NODES[e].id.begin(), NODES[e].id.end(), 'r') == NODES[e].id.end() || NODES[e].id[0] == 's') {
			out3 << NODES[e].id << ",";
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					out3 << allReactiveLoad[i][e][k] << ",";
				}
			}
			out3 << endl;
		}
	}
	out3.close();
	output_file = gridlabd_path + string("/generatorsReal_") + ss.str() + string(".txt");
	ofstream out4(output_file.c_str());
	for (int j = 0; j < NumNodes; ++j) {
		int e = voltageSum[j].ind;
		if (find(NODES[e].id.begin(), NODES[e].id.end(), 'r') == NODES[e].id.end() || NODES[e].id[0] == 's') {
			out4 << NODES[e].id << ",";
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					out4 << allRealGen[i][e][k] << ",";
				}
			}
			out4 << endl;
		}
	}
	out4.close();
	output_file = gridlabd_path + string("/generatorsReactive_") + ss.str() + string(".txt");
	ofstream out5(output_file.c_str());
	for (int j = 0; j < NumNodes; ++j) {
		int e = voltageSum[j].ind;
		if (find(NODES[e].id.begin(), NODES[e].id.end(), 'r') == NODES[e].id.end() || NODES[e].id[0] == 's') {
			out5 << NODES[e].id << ",";
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					out5 << allReactiveGen[i][e][k] << ",";
				}
			}
			out5 << endl;
		}
	}
	out5.close();
	output_file = gridlabd_path + string("/lines_") + ss.str() + string(".txt");
	ofstream out6(output_file.c_str());
	for (int e = 0; e < NumEdges; ++e) {
		out6 << EDGES[e].id << "," << EDGES[e].node1id << "," << EDGES[e].node2id
		     << "," << EDGES[e].length << "," << EDGES[e].capacity;
		for (int k = 0; k < NumPhases; ++k) {
			if (EDGES[e].hasphase[k]) {
				out6 << ",1";
			} else {
				out6 << ",0";
			}
		}
		if (!EDGES[e].istransformer) {
			map<int, lineCodeData>::iterator it = LINECODES.find(EDGES[e].linecode);

			out6 << "," << it->first;
			for (int k = 0; k < EDGES[e].NumPhases; ++k) {
				for (int l = 0; l < EDGES[e].NumPhases; ++l) {
					out6 << "," << (it->second.rmatrix[0][k][l]);
				}

			}
			for (int k = 0; k < EDGES[e].NumPhases; ++k) {
				for (int l = 0; l < EDGES[e].NumPhases; ++l) {
					out6 << "," << (it->second.xmatrix[0][k][l]);
				}

			}
		}
		out6 << endl;
	}
	out6.close();

}


// FOLLOWING METHOD IS USED TO SOLVE BENDERS DECOMPOSITION
void BendersFullScenarioGeneration(
	vector<double> &ObjectiveVector,
	double MIN_FEAS_OBJ,
	vector<double> &MIN_FEAS_SOL,
	vector<IloModel> &SubProblem,
	vector<Solver> &SubSolver,
	map<string, int> &hashTableVertex,
	map<string, int> &hashTableEdge,
	map<string, int> &hashTableGenerator,
	map<string, int> &hashTableLoad,
	vector<nodeData> &NODES,
	vector<edgeData> &EDGES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	map<int, lineCodeData> &LINECODES,
	vector<dataNode<double> > &LINE_CONSTRUCTION_COST,
	vector<dataNode<double> > &MICROGRID_COST,
	vector<dataNode<double> > &MICROGRID_FIXED_COST,
	vector<dataNode<bool> > &IS_CRITICAL_LOAD,
	vector<dataNode<double> > &MAX_MICROGRID,
	vector<dataNode<double> > &HARDEN_COST,
	vector<dataNode<double> > &LINE_SWITCH_COST,
	double LoadMet,
	double CriticalLoadMet,
	double PhaseVariation,
	vector<graph_vertex> &G,
	vector<vector<int> > &CYCLES,
	vector<IloObjective> &SubObjective, // OBJECTIVE 
	vector<IloRangeArray> &lineExistsConstraint, // CONSTRAINTS
	vector<IloRangeArray> &lineConstraint1, // CONSTRAINTS
	vector<IloRangeArray> &lineConstraint2,
	vector<IloRangeArray> &lineReactiveConstraint1, // CONSTRAINTS
	vector<IloRangeArray> &lineReactiveConstraint2,
	vector<IloRangeArray> &directionConstraint,
	vector<IloRangeArray> &switchConstraint1,
	vector<IloRangeArray> &switchConstraint2,
	vector<IloRangeArray> &switchReactiveConstraint1,
	vector<IloRangeArray> &switchReactiveConstraint2,
	vector<IloRangeArray> &variationConstraint1,
	vector<IloRangeArray> &variationConstraint2,
	vector<IloRangeArray> &variationReactiveConstraint1,
	vector<IloRangeArray> &variationReactiveConstraint2,
	vector<IloRangeArray> &damageConstraint,
	vector<IloRangeArray> &loadConstraint,
	vector<IloRangeArray> &loadReactiveConstraint,
	vector<IloRangeArray> &generationConstraint,
	vector<IloRangeArray> &generationReactiveConstraint,
	vector<IloRangeArray> &balanceConstraint,
	vector<IloRangeArray> &balanceReactiveConstraint,
	vector<IloRangeArray> &microgridConstraint,
	vector<IloRangeArray> &cycleConstraint,
	vector<IloRangeArray> &cycleSubConstraint,
	vector<IloRangeArray> &switchableConstraint,
	vector<IloRangeArray> &criticalServeConstraint,
	vector<IloRangeArray> &criticalServeReactiveConstraint,
	vector<IloRangeArray> &loadServeConstraint,
	vector<IloRangeArray> &loadServeReactiveConstraint,
	vector<IloRangeArray> &switchImplicationConstraint1,
	//vector<IloRangeArray> &switchImplicationConstraint2,
	vector<IloRangeArray> &hardenConstraint,
	vector<IloRangeArray> &voltageConstraint,
	vector<IloRangeArray> &voltageOffsetConstraint,
	vector<IloRangeArray> &voltageVariationConstraint1,
	vector<IloRangeArray> &voltageVariationConstraint2,
	vector<IloNumVarArray> &lineUseVariable, // VARIABLES
	vector<IloNumVarArray> &lineExistsVariable, // VARIABLES
	vector<IloNumVarArray> &lineDirectionVariable,
	vector<IloNumVarArray> &lineCycleVariable,
	vector<IloNumVarArray> &lineHardenVariable,
	vector<IloNumVarArray> &switchUseVariable,
	vector<IloNumVarArray> &switchCycleVariable,
	vector<IloNumVarArray> &flowVariable,
	vector<IloNumVarArray> &flowReactiveVariable,
	vector<IloNumVarArray> &voltageVariable,
	vector<IloNumVarArray> &voltageOffsetVariable,
	vector<IloNumVarArray> &loadVariable,
	vector<IloNumVarArray> &loadReactiveVariable,
	vector<IloNumVarArray> &loadServeVariable,
	vector<IloNumVarArray> &generatorVariable,
	vector<IloNumVarArray> &generatorReactiveVariable,
	vector<IloNumVarArray> &facilityVariable,
	const vector<int> &NumDamagedEdges,
	int NumPhases,
	int NumEdges,
	int NumLoads,
	int NumGenerators,
	int NumScenarios,
	double discreteMicrogrid,
	string damage_str,
	string extent_str,
	string root_str,
	bool HardenedDisabled,
	bool ExactSolver,
	bool LDFIndicator,
	int MaxVoltageRounds,
	bool CHECK_RESULTS,
	vector<vector<dataNode<bool> > > &DISABLED,
	vector<vector<dataNode<bool> > > &HARDENED_DISABLED,
	double newImpedanceMultiplier) {
	
	bool isInfeasible = true;	
	int numCon = -1;
	int NumNodes = NODES.size();

	vector<double> IterationObjective(0);
	vector<double> IterationCPUTIME(0);

	// SET THE OBJECTIVE VALUE OF SUB PROBLEMS:
	vector<vector<double> > TotalDemandPerPhase(NumScenarios, vector<double>(NumPhases, 0.0));
	vector<IloNumVar> minInfeasibilityVariable(NumScenarios);
	for (int i = 0; i < NumScenarios; ++i) {
		IloEnv subEnv = SubSolver[i].cplex.getEnv();
		SubProblem[i].remove(SubObjective[i]);
		SubObjective[i].end();
		SubObjective[i] = IloObjective(subEnv);
		SubObjective[i].setSense(IloObjective::Maximize);


		string str = "minInfeasibilityVariable";
		minInfeasibilityVariable[i] = IloNumVar(subEnv, 0.0, IloInfinity, ILOFLOAT, str.c_str());
		criticalServeConstraint[i][0].setLinearCoef(minInfeasibilityVariable[i], 1.0);
		criticalServeReactiveConstraint[i][0].setLinearCoef(minInfeasibilityVariable[i], 1.0);
		loadServeConstraint[i][0].setLinearCoef(minInfeasibilityVariable[i], 1.0);
		loadServeReactiveConstraint[i][0].setLinearCoef(minInfeasibilityVariable[i], 1.0);
		SubObjective[i].setLinearCoef(minInfeasibilityVariable[i], -1.0);
		SubProblem[i].add(minInfeasibilityVariable[i]);
				
		/*SubProblem[i].remove(criticalServeConstraint[i]);
		for (int j = 0; j < IS_CRITICAL_LOAD.size(); ++j) {
			if (IS_CRITICAL_LOAD[j].data) {
				for (int k = 0; k < NumPhases; ++k) {
					if (LOADS[hashTableLoad[IS_CRITICAL_LOAD[j].id]].hasphase[k]) {
						TotalDemandPerPhase[i][k] += LOADS[hashTableLoad[IS_CRITICAL_LOAD[j].id]].realphase[k];
						SubObjective[i].setLinearCoef(loadVariable[i][hashTableVertex[LOADS[hashTableLoad[IS_CRITICAL_LOAD[j].id]].node_id] * NumPhases + k], 1.0);
					}
				}
			}
		}*/
		
		SubProblem[i].add(SubObjective[i]);
		//SubSolver[i].cplex.exportModel("subproblem.lp");
		
		// REMOVE FACILITIES FROM SUBPROBLEM
		SubProblem[i].remove(facilityVariable[i]);
		SubProblem[i].remove(generationConstraint[i]);
		SubProblem[i].remove(generationReactiveConstraint[i]);
			
	}

	int ColSize = 3 * NumEdges + NumGenerators;

	vector<vector<vector<double> > > subInfeasibilityCuts(NumScenarios, vector<vector<double> >(0));

	IloEnv env = SubSolver[NumScenarios].cplex.getEnv();

	vector<vector<double> > InfCutCoef;
	vector<double> InfCutRHS;

	IloRangeArray InfeasibilityCuts(env);
	int numInfeasibilityCuts = -1;

	IloRangeArray BaseLineVarConstraint(env);
	IloRangeArray BaseSwitchVarConstraint(env);
	IloRangeArray BaseHardenVarConstraint(env);
	IloRangeArray BaseMicrogridVarConstraint(env);
	IloRangeArray BaseFacilityVarConstraint(env);
	IloRangeArray BaseFacilityVarReactiveConstraint(env);

	SubProblem[NumScenarios].add(InfeasibilityCuts);
	SubProblem[NumScenarios].add(BaseLineVarConstraint);
	SubProblem[NumScenarios].add(BaseSwitchVarConstraint);
	SubProblem[NumScenarios].add(BaseHardenVarConstraint);
	SubProblem[NumScenarios].add(BaseMicrogridVarConstraint);
	SubProblem[NumScenarios].add(BaseFacilityVarConstraint);
	SubProblem[NumScenarios].add(BaseFacilityVarReactiveConstraint);

	SubProblem[NumScenarios].remove(lineExistsConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(lineConstraint1[NumScenarios]);
	SubProblem[NumScenarios].remove(lineConstraint2[NumScenarios]);
	SubProblem[NumScenarios].remove(lineReactiveConstraint1[NumScenarios]);
	SubProblem[NumScenarios].remove(lineReactiveConstraint2[NumScenarios]);
	SubProblem[NumScenarios].remove(directionConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(switchConstraint1[NumScenarios]);
	SubProblem[NumScenarios].remove(switchConstraint2[NumScenarios]);
	SubProblem[NumScenarios].remove(switchReactiveConstraint1[NumScenarios]);
	SubProblem[NumScenarios].remove(switchReactiveConstraint2[NumScenarios]);
	SubProblem[NumScenarios].remove(variationConstraint1[NumScenarios]);
	SubProblem[NumScenarios].remove(variationConstraint2[NumScenarios]);
	SubProblem[NumScenarios].remove(variationReactiveConstraint1[NumScenarios]);
	SubProblem[NumScenarios].remove(variationReactiveConstraint2[NumScenarios]);
	SubProblem[NumScenarios].remove(damageConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(loadConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(loadReactiveConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(generationConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(generationReactiveConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(balanceConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(balanceReactiveConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(microgridConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(cycleConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(cycleSubConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(switchableConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(criticalServeConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(criticalServeReactiveConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(loadServeConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(loadServeReactiveConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(switchImplicationConstraint1[NumScenarios]);
	SubProblem[NumScenarios].remove(hardenConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(voltageConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(voltageOffsetConstraint[NumScenarios]);
	SubProblem[NumScenarios].remove(voltageVariationConstraint1[NumScenarios]);
	SubProblem[NumScenarios].remove(voltageVariationConstraint2[NumScenarios]);

	IloObjective cloneSubObjective = IloObjective(env);

	vector<IloRangeArray> clonelineExistsConstraint(0);
	vector<IloRangeArray> clonelineConstraint1(0);
	vector<IloRangeArray> clonelineConstraint2(0);
	vector<IloRangeArray> clonelineReactiveConstraint1(0);
	vector<IloRangeArray> clonelineReactiveConstraint2(0);
	vector<IloRangeArray> clonedirectionConstraint(0);
	vector<IloRangeArray> cloneswitchConstraint1(0);
	vector<IloRangeArray> cloneswitchConstraint2(0);
	vector<IloRangeArray> cloneswitchReactiveConstraint1(0);
	vector<IloRangeArray> cloneswitchReactiveConstraint2(0);
	vector<IloRangeArray> clonevariationConstraint1(0);
	vector<IloRangeArray> clonevariationConstraint2(0);
	vector<IloRangeArray> clonevariationReactiveConstraint1(0);
	vector<IloRangeArray> clonevariationReactiveConstraint2(0);
	vector<IloRangeArray> clonedamageConstraint(0);
	vector<IloRangeArray> cloneloadConstraint(0);
	vector<IloRangeArray> cloneloadReactiveConstraint(0);
	vector<IloRangeArray> clonegenerationConstraint(0);
	vector<IloRangeArray> clonegenerationReactiveConstraint(0);
	vector<IloRangeArray> clonebalanceConstraint(0);
	vector<IloRangeArray> clonebalanceReactiveConstraint(0);
	vector<IloRangeArray> clonemicrogridConstraint(0);
	vector<IloRangeArray> clonecycleConstraint(0);
	vector<IloRangeArray> clonecycleSubConstraint(0);
	vector<IloRangeArray> cloneswitchableConstraint(0);
	vector<IloRangeArray> clonecriticalServeConstraint(0);
	vector<IloRangeArray> clonecriticalServeReactiveConstraint(0);
	vector<IloRangeArray> cloneloadServeConstraint(0);
	vector<IloRangeArray> cloneloadServeReactiveConstraint(0);
	vector<IloRangeArray> cloneswitchImplicationConstraint1(0);
	//vector<IloRangeArray> cloneswitchImplicationConstraint2(0);
	vector<IloRangeArray> clonehardenConstraint(0);
	vector<IloRangeArray> clonevoltageConstraint(0);
	vector<IloRangeArray> clonevoltageOffsetConstraint(0);
	vector<IloRangeArray> clonevoltageVariationConstraint1(0);
	vector<IloRangeArray> clonevoltageVariationConstraint2(0);

	// VOLTAGE INFEASIBILITY CONSTRAINTS
	// vector<IloRangeArray> clonevoltageInfeasibilityConstraint(0);
	// VOLTAGE INFEASIBILITY CUT DATA
	vector<int> scenarioList;
	vector<list<string> > nodesList;
	vector<string> infNodeList;
	vector<list<string> > edgesList;


	vector<IloNumVarArray> clonelineUseVariable(0);
	vector<IloNumVarArray> clonelineExistsVariable(0);
	vector<IloNumVarArray> clonelineDirectionVariable(0);
	vector<IloNumVarArray> clonelineCycleVariable(0);
	vector<IloNumVarArray> clonelineHardenVariable(0);
	vector<IloNumVarArray> cloneswitchUseVariable(0);
	vector<IloNumVarArray> cloneswitchCycleVariable(0);
	vector<IloNumVarArray> cloneflowVariable(0);
	vector<IloNumVarArray> cloneflowReactiveVariable(0);
	vector<IloNumVarArray> clonevoltageVariable(0);
	vector<IloNumVarArray> clonevoltageOffsetVariable(0);
	vector<IloNumVarArray> cloneloadVariable(0);
	vector<IloNumVarArray> cloneloadReactiveVariable(0);
	vector<IloNumVarArray> cloneloadServeVariable(0);
	vector<IloNumVarArray> clonegeneratorVariable(0);
	vector<IloNumVarArray> clonegeneratorReactiveVariable(0);
	vector<IloNumVarArray> clonefacilityVariable(0);

	vector<IloNumArray> cloneFlowVals(0);
	vector<IloNumArray> cloneFlowReactiveVals(0);

	vector<IloNumArray> cloneLineExistsVals(0);
	vector<IloNumArray> cloneLineUseVals(0);
	vector<IloNumArray> cloneLineDirectionVals(0);
	vector<IloNumArray> cloneSwitchUseVals(0);
	vector<IloNumArray> cloneLineHardenVals(0);
	vector<IloNumArray> cloneGeneratorVals(0);
	vector<IloNumArray> cloneGeneratorReactiveVals(0);
	vector<IloNumArray> cloneLoadVals(0);
	vector<IloNumArray> cloneLoadReactiveVals(0);
	vector<IloNumArray> cloneLoadServeVals(0);
	vector<IloNumArray> cloneVoltageVals(0);


	// THESE VECTORS WILL GATHER SOLUTION STATS PER SCENARIO
	vector<IloNumArray> subLineUseVals(NumScenarios);
	vector<IloNumArray> subLineDirectionVals(NumScenarios);
	vector<IloNumArray> subSwitchUseVals(NumScenarios);
	vector<IloNumArray> subGeneratorVals(NumScenarios);
	vector<IloNumArray> subGeneratorReactiveVals(NumScenarios);
	vector<IloNumArray> subLoadVals(NumScenarios);
	vector<IloNumArray> subLoadReactiveVals(NumScenarios);
	vector<IloNumArray> subLoadServeVals(NumScenarios);
	vector<IloNumArray> subVoltageVals(NumScenarios);
	for (int i = 0; i < NumScenarios; ++i) {
		IloEnv subEnv = SubSolver[i].cplex.getEnv();
		subLineUseVals[i] = IloNumArray(subEnv);
		subLineDirectionVals[i] = IloNumArray(subEnv);
		subSwitchUseVals[i] = IloNumArray(subEnv);
		subGeneratorVals[i] = IloNumArray(subEnv);
		subGeneratorReactiveVals[i] = IloNumArray(subEnv);
		subLoadVals[i] = IloNumArray(subEnv);
		subLoadReactiveVals[i] = IloNumArray(subEnv);
		subLoadServeVals[i] = IloNumArray(subEnv);
		subVoltageVals[i] = IloNumArray(subEnv);
	}


	cout << "NUM DAMAGED STATS:" << endl;
	cout << "Max element = " << maxelement<int>(NumDamagedEdges) << endl;
	cout << "Min element = " << minelement<int>(NumDamagedEdges) << endl;
	/*for (int i = 0; i < NumDamagedEdges.size(); ++i) {
		cout << NumDamagedEdges[i] << " ";
	}
	cout << endl;*/


	vector<vector<double> > MAXREALPHASE(NODES.size(), vector<double>(3, 0.0));
	vector<vector<double> > MAXREACTIVEPHASE(NODES.size(), vector<double>(3, 0.0));
	for (int j = 0; j < NODES.size(); ++j) {
		for (int k = 0; k < GENERATORS.size(); ++k) {
			for (int l = 0; l < NumPhases; ++l) {
				if (hashTableVertex[GENERATORS[k].node_id] == j && GENERATORS[k].hasphase[l]) {
					MAXREALPHASE[j][l] += GENERATORS[k].maxrealphase[l];
					MAXREACTIVEPHASE[j][l] += GENERATORS[k].maxreactivephase[l];
				}
			}
		}
	}
	vector<double> MAX_ADDED(GENERATORS.size(), 0.0);
	for (int j = 0; j < MAX_MICROGRID.size(); ++j) {
		MAX_ADDED[hashTableGenerator[MAX_MICROGRID[j].id]] = MAX_MICROGRID[j].data;
	}

	SubSolver[NumScenarios].cplex.setParam(IloCplex::ClockType, 1);

	double MAX_CPUTIME = 72.0 * 60.0 * 60.0;
	bool FirstVoltageCut = true;
	double CPUTIME_SBD = MAX_CPUTIME;

	vector<bool> voltageInfeasibleScenarios(NumScenarios, false);
	int bendersIter = -1;
	int MAX_SCENARIOS_ADDED = NumScenarios;//NumScenarios;
	bool ADD_INFEASIBILITY_CUTS = false;
	vector<bool> ScenarioIncluded(NumScenarios, false);
	double incumbentObjective;
	time_t begin = clock();
	bool TIMED_WITHOUT_SOLUTION = false;
	string STATUS = "Feasible";
	vector<double> incumbentSolution(ColSize, 0.0);
	vector<double> prevSolution(ColSize, 0.0);
	vector<int> addedScenarios(0);
	set<vector<int> > optimalSolutions;
	vector<set<vector<int> > > componentSolutions(NumScenarios, set<vector<int> >());
	int numVoltageCuts = 0;
	do {
		IloNumArray lineusevals(env);
		IloNumArray lineswitchvals(env);
		IloNumArray linehardenvals(env);
		IloNumArray facilityvals(env);

		// TODO: Implement a local search method here.
		for (int j = 0; j < ColSize; ++j) {
			incumbentSolution[j] = MIN_FEAS_SOL[j];
		}
		incumbentObjective = std::inner_product(ObjectiveVector.begin(), ObjectiveVector.end(), incumbentSolution.begin(), 0.0);

		cout << "CPUTIME: " << ((double) (clock() - begin) / (double) CLOCKS_PER_SEC) << endl;
		if (ExactSolver) {
			SubSolver[NumScenarios].cplex.setParam(IloCplex::TiLim, MAX_CPUTIME - ((double) (clock() - begin) / (double) CLOCKS_PER_SEC) > 60.0 ? MAX_CPUTIME - ((double) (clock() - begin) / (double) CLOCKS_PER_SEC) : 60.0);
			SubSolver[NumScenarios].cplex.solve();
			//cout << "SOLVED " << SubSolver[NumScenarios].cplex.getStatus() << endl;

			//cout << "CHECK SOLUTION" << endl;
			if (SubSolver[NumScenarios].cplex.getStatus() == IloAlgorithm::Feasible || SubSolver[NumScenarios].cplex.getStatus() == IloAlgorithm::Optimal) {
				SubSolver[NumScenarios].cplex.getValues(lineusevals, lineUseVariable[NumScenarios]);
				SubSolver[NumScenarios].cplex.getValues(lineswitchvals, switchUseVariable[NumScenarios]);
				SubSolver[NumScenarios].cplex.getValues(linehardenvals, lineHardenVariable[NumScenarios]);
				SubSolver[NumScenarios].cplex.getValues(facilityvals, facilityVariable[NumScenarios]);
				//incumbentObjective = SubSolver[NumScenarios].cplex.getObjValue();
				
				for (int i = 0; i < NumEdges; ++i) {
					incumbentSolution[i] = lineusevals[i];
					incumbentSolution[i+NumEdges] = lineswitchvals[i];
					incumbentSolution[i+2*NumEdges] = linehardenvals[i];
				}
				for (int j = 0; j < NumGenerators; ++j) {
					incumbentSolution[3*NumEdges + j] = facilityvals[j];
				}
				incumbentObjective = std::inner_product(ObjectiveVector.begin(), ObjectiveVector.end(), incumbentSolution.begin(), 0.0);
				STATUS = "Optimal";

				cloneFlowVals.clear();
				cloneFlowReactiveVals.clear();

				cloneLineExistsVals.clear();
				cloneLineUseVals.clear();
				cloneLineDirectionVals.clear();
				cloneSwitchUseVals.clear();
				cloneLineHardenVals.clear();
				cloneGeneratorVals.clear();
				cloneGeneratorReactiveVals.clear();
				cloneLoadVals.clear();
				cloneLoadReactiveVals.clear();
				cloneLoadServeVals.clear();
				cloneVoltageVals.clear();
				for (int j = 0; j < (bendersIter+1); ++j) {
					cloneFlowVals.push_back(IloNumArray(env));
					cloneFlowReactiveVals.push_back(IloNumArray(env));
					cloneLineExistsVals.push_back(IloNumArray(env));
					cloneLineUseVals.push_back(IloNumArray(env));
					cloneLineHardenVals.push_back(IloNumArray(env));
					cloneLineDirectionVals.push_back(IloNumArray(env));
					cloneSwitchUseVals.push_back(IloNumArray(env));
					cloneGeneratorVals.push_back(IloNumArray(env));
					cloneGeneratorReactiveVals.push_back(IloNumArray(env));
					cloneLoadVals.push_back(IloNumArray(env));
					cloneLoadReactiveVals.push_back(IloNumArray(env));
					cloneLoadServeVals.push_back(IloNumArray(env));
					cloneVoltageVals.push_back(IloNumArray(env));
					//cout << clonelineUseVariable[j].getSize() << endl;
					SubSolver[NumScenarios].cplex.getValues(cloneFlowVals[j], cloneflowVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneFlowReactiveVals[j], cloneflowReactiveVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLineExistsVals[j], clonelineExistsVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLineUseVals[j], clonelineUseVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLineHardenVals[j], clonelineHardenVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLineDirectionVals[j], clonelineDirectionVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneSwitchUseVals[j], cloneswitchUseVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneGeneratorVals[j], clonegeneratorVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneGeneratorReactiveVals[j], clonegeneratorReactiveVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLoadVals[j], cloneloadVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLoadReactiveVals[j], cloneloadReactiveVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneLoadServeVals[j], cloneloadServeVariable[j]);
					SubSolver[NumScenarios].cplex.getValues(cloneVoltageVals[j], clonevoltageVariable[j]);
				}
	
			} else {
				cout << (SubSolver[NumScenarios].cplex.getStatus() == IloAlgorithm::Infeasible ? "INFEASIBLE PROBLEM" : "UNBOUNDED PROBLEM") <<  endl;
				TIMED_WITHOUT_SOLUTION = true;
				STATUS = "Feasible";
				break;
			}
		} else {
			// -- GET LP RELAXATION ----------------------------------------------------------------
			vector<double> LPSolution(ColSize, 0.0);

			solveRelaxed(
				SubProblem[NumScenarios],
				SubSolver[NumScenarios],
				lineUseVariable[NumScenarios], // VARIABLES
				lineExistsVariable[NumScenarios], // VARIABLES
				lineDirectionVariable[NumScenarios],
				lineCycleVariable[NumScenarios],
				lineHardenVariable[NumScenarios],
				switchUseVariable[NumScenarios],
				switchCycleVariable[NumScenarios],
				flowVariable[NumScenarios],
				flowReactiveVariable[NumScenarios],
				voltageVariable[NumScenarios],
				voltageOffsetVariable[NumScenarios],
				loadVariable[NumScenarios],
				loadReactiveVariable[NumScenarios],
				loadServeVariable[NumScenarios],
				generatorVariable[NumScenarios],
				generatorReactiveVariable[NumScenarios],
				facilityVariable[NumScenarios],
				clonelineUseVariable, // SUB PROBLEM VARIABLES
				clonelineExistsVariable, // SUB PROBLEM VARIABLES
				clonelineDirectionVariable,
				clonelineCycleVariable,
				clonelineHardenVariable,
				cloneswitchUseVariable,
				cloneswitchCycleVariable,
				cloneflowVariable,
				cloneflowReactiveVariable,
				clonevoltageVariable,
				clonevoltageOffsetVariable,
				cloneloadVariable,
				cloneloadReactiveVariable,
				cloneloadServeVariable,
				clonegeneratorVariable,
				clonegeneratorReactiveVariable,
				LPSolution,
				NumScenarios,
				NumEdges,
				NumGenerators,
				NumPhases);
			// -------------------------------------------------------------------------------------
			

			VariableNeighborhoodSearchSolver(
				MIN_FEAS_OBJ,
				incumbentSolution,
				LPSolution,
				SubProblem,
				SubSolver,
				hashTableVertex,
				hashTableEdge,
				hashTableGenerator,
				hashTableLoad,
				NODES,
				EDGES,
				GENERATORS,
				LOADS,
				LINECODES,
				LINE_CONSTRUCTION_COST,
				MICROGRID_COST,
				MICROGRID_FIXED_COST,
				IS_CRITICAL_LOAD,
				MAX_MICROGRID,
				HARDEN_COST,
				LINE_SWITCH_COST,
				LoadMet,
				CriticalLoadMet,
				PhaseVariation,
				G,
				CYCLES,
				lineUseVariable, // VARIABLES
				lineExistsVariable, // VARIABLES
				lineDirectionVariable,
				lineCycleVariable,
				lineHardenVariable,
				switchUseVariable,
				switchCycleVariable,
				flowVariable,
				flowReactiveVariable,
				voltageVariable,
				voltageOffsetVariable,
				loadVariable,
				loadReactiveVariable,
				loadServeVariable,
				generatorVariable,
				generatorReactiveVariable,
				facilityVariable,
				clonelineUseVariable,
				clonelineExistsVariable,
				clonelineDirectionVariable,
				cloneswitchUseVariable,
				clonegeneratorVariable,
				clonegeneratorReactiveVariable,
				clonevoltageVariable,
				clonevoltageOffsetVariable,
				cloneloadVariable,
				cloneloadReactiveVariable,
				cloneloadServeVariable,
				cloneLineExistsVals,
				cloneLineUseVals,
				cloneLineDirectionVals,
				cloneSwitchUseVals,
				cloneGeneratorVals,
				cloneGeneratorReactiveVals,
				cloneLoadVals,
				cloneLoadReactiveVals,
				cloneLoadServeVals,
				cloneVoltageVals,
				bendersIter,
				NumDamagedEdges,
				NumPhases,
				NumEdges,
				NumGenerators,
				NumScenarios,
				discreteMicrogrid,
				damage_str,
				extent_str,
				root_str,
				HardenedDisabled,
				MAX_CPUTIME - ((double) (clock() - begin) / (double) CLOCKS_PER_SEC) > 60.0 ? MAX_CPUTIME - ((double) (clock() - begin) / (double) CLOCKS_PER_SEC) : 60.0);

			STATUS = "Feasible";
			
			for (int j = 0; j < NumEdges; ++j) {
				lineusevals.add(incumbentSolution[j]);
				lineswitchvals.add(incumbentSolution[j + NumEdges]);
				linehardenvals.add(incumbentSolution[j + 2*NumEdges]);
			}
			for (int j = 0; j < NumGenerators; ++j) {
				facilityvals.add(incumbentSolution[j + 3*NumEdges]);
			}

			incumbentObjective = std::inner_product(ObjectiveVector.begin(), ObjectiveVector.end(), incumbentSolution.begin(), 0.0);
	
		}


		IterationObjective.push_back(incumbentObjective);
		IterationCPUTIME.push_back(((double) (clock() - begin) / (double) CLOCKS_PER_SEC));

		vector<sortWithIndices> CurrentObjective(NumScenarios);
		
		isInfeasible = false;

		vector<bool> isSubInfeasible(NumScenarios, false);

		#pragma omp parallel for num_threads(32)
		for (int i = 0; i < NumScenarios; ++i) {
			if (!ScenarioIncluded[i]) {
				//cout << "THRED #" << omp_get_thread_num() << endl;
				CurrentObjective[i].ind = i;
				CurrentObjective[i].tie = NumDamagedEdges[i];

				vector<double> lineUseUB(NumEdges);
				vector<double> lineUseLB(NumEdges);
				vector<double> switchUseUB(NumEdges);
				vector<double> lineHardenUB(NumEdges);
				vector<double> lineHardenLB(NumEdges);
				vector<double> generatorUB(NumNodes*NumPhases);
				vector<double> generatorUBReactive(NumNodes*NumPhases);
				// MODIFY THE ORIGINAL SUBPROBLEM:
				for (int j = 0; j < NumEdges; ++j) {
					lineUseUB[j] = lineUseVariable[i][j].getUB();
					lineUseLB[j] = lineUseVariable[i][j].getLB();
					switchUseUB[j] = switchUseVariable[i][j].getUB();
					lineHardenUB[j] = lineHardenVariable[i][j].getUB();
					lineHardenLB[j] = lineHardenVariable[i][j].getLB();
					lineUseVariable[i][j].setUB(round((double)lineusevals[j]));
					if (DISABLED[i][j].data == false) {
						lineUseVariable[i][j].setLB(round((double)lineusevals[j]));
					}
					switchUseVariable[i][j].setUB(round((double)lineswitchvals[j]));
					// should be checked if HARDENED_DISABLED is implemented as an upper bound instead of constraint
					lineHardenVariable[i][j].setUB(lineHardenUB[j] < round((double)linehardenvals[j]) ? lineHardenUB[j] : round((double)linehardenvals[j]));
					if (HARDENED_DISABLED[i][j].data == false) {
						lineHardenVariable[i][j].setLB(lineHardenUB[j] < round((double)linehardenvals[j]) ? lineHardenUB[j] : round((double)linehardenvals[j]));
					}
				}
				for (int j = 0; j < NumNodes; ++j) {
					for (int k = 0; k < NumPhases; ++k) {
						generatorUB[j * NumPhases + k] = generatorVariable[i][j * NumPhases + k].getUB();
						generatorUBReactive[j * NumPhases + k] = generatorReactiveVariable[i][j * NumPhases + k].getUB();
						double newUB = MAXREALPHASE[j][k];
						double newUBReactive = MAXREACTIVEPHASE[j][k];
						for (int l = 0; l < NumGenerators; ++l) {
							if (hashTableVertex[GENERATORS[l].node_id] == j && GENERATORS[l].hasphase[k]) {
								newUB += MAX_ADDED[l]*round((double)facilityvals[l]);	
								newUBReactive += MAX_ADDED[l]*round((double)facilityvals[l]);	
							}
						}
						generatorVariable[i][j * NumPhases + k].setUB(newUB > 0.0 ? newUB : 0.0);
						generatorReactiveVariable[i][j * NumPhases + k].setUB(newUBReactive > 0.0 ? newUBReactive : 0.0);
					}
				}
				SubSolver[i].cplex.solve();

				if (SubSolver[i].cplex.getStatus() == IloAlgorithm::Infeasible || SubSolver[i].cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded) {
					//CurrentObjective[i].val = 0.0 - CriticalLoadMet * vectorsum(TotalDemandPerPhase[i]);
					CurrentObjective[i].val = -1e+5;
				} else {
					CurrentObjective[i].val = SubSolver[i].cplex.getObjValue();
					//CurrentObjective[i].val = SubSolver[i].cplex.getObjValue() - CriticalLoadMet * vectorsum(TotalDemandPerPhase[i]); 
					SubSolver[i].cplex.getValues(subLineUseVals[i], lineUseVariable[i]);
					SubSolver[i].cplex.getValues(subLineDirectionVals[i], lineDirectionVariable[i]);
					SubSolver[i].cplex.getValues(subSwitchUseVals[i], switchUseVariable[i]);
					SubSolver[i].cplex.getValues(subGeneratorVals[i], generatorVariable[i]);
					SubSolver[i].cplex.getValues(subGeneratorReactiveVals[i], generatorReactiveVariable[i]);
					SubSolver[i].cplex.getValues(subLoadVals[i], loadVariable[i]);
					SubSolver[i].cplex.getValues(subLoadReactiveVals[i], loadReactiveVariable[i]);
					SubSolver[i].cplex.getValues(subLoadServeVals[i], loadServeVariable[i]);
					SubSolver[i].cplex.getValues(subVoltageVals[i], voltageVariable[i]);
				}	
				// RECOVER THE ORIGINAL SUBPROBLEM:
				for (int j = 0; j < NumEdges; ++j) {
					lineUseVariable[i][j].setUB(lineUseUB[j]);
					lineUseVariable[i][j].setLB(lineUseLB[j]);
					switchUseVariable[i][j].setUB(switchUseUB[j]);
					lineHardenVariable[i][j].setUB(lineHardenUB[j]);
					lineHardenVariable[i][j].setLB(lineHardenLB[j]);
				}
				for (int j = 0; j < NumNodes; ++j) {
					for (int k = 0; k < NumPhases; ++k) {
						generatorVariable[i][j * NumPhases + k].setUB(generatorUB[j * NumPhases + k]);
						generatorReactiveVariable[i][j * NumPhases + k].setUB(generatorUBReactive[j * NumPhases + k]);
					}
				}
			}
		}
		std::sort(CurrentObjective.begin(), CurrentObjective.end(), mycomperator);
		cout << "Current most infeasible objective: " << CurrentObjective[0].val << endl;

		cout << "CURRENT STATUS: " << SubSolver[NumScenarios].cplex.getStatus() << endl;

		int numEdgesUsed = 0;
		int numSwitchesUsed = 0;
		int numHardenUsed = 0;
		int numFacilityUsed = 0;
		//cout << "LINES: ";
		for (int j = 0; j < NumEdges; ++j) {
			for (int l = 0; l < LINE_CONSTRUCTION_COST.size(); ++l) {
				if (fabs((double)lineusevals[j] - 1.0) < 1e-4 && hashTableEdge[LINE_CONSTRUCTION_COST[l].id] == j && LINE_CONSTRUCTION_COST[l].data > 0.0) {
					numEdgesUsed++;
					break;
					//cout << EDGES[j].id << " ";
				}
			}
		}
		//cout << endl;
		//cout << "SWITCHES: ";
		for (int j = 0; j < NumEdges; ++j) {
			for (int l = 0; l < LINE_SWITCH_COST.size(); ++l) {
				if (fabs((double)lineswitchvals[j] - 1.0) < 1e-4  && hashTableEdge[LINE_SWITCH_COST[l].id] == j && LINE_SWITCH_COST[l].data > 0.0 ) {
					numSwitchesUsed++;
					//cout << EDGES[j].id << " ";
				}
			}
		}
		//cout << endl;
		//cout << "HARDENED: ";
		for (int j = 0; j < NumEdges; ++j) {
			if (fabs((double)linehardenvals[j] - 1.0) < 1e-4) {
				numHardenUsed++;
				//cout << EDGES[j].id << " ";
			}
		}
		//cout << endl;

		/*vector<int> incumbentIntSolution(ColSize);
		for (int j = 0; j < ColSize; ++j) {
			incumbentIntSolution[j] = (int) round(incumbentSolution[j]);
		}
		if (optimalSolutions.find(incumbentIntSolution) == optimalSolutions.end()) optimalSolutions.insert(incumbentIntSolution);
		else {cout << "********* PREVIOUS SOLUTION **********" << endl; return;}*/

		for (int j = 0; j < NumGenerators; ++j) {
			if (fabs((double)facilityvals[j] - 1.0) < 1e-4) numFacilityUsed++;
		}
		cout << "NUM NEW EDGES: " << numEdgesUsed << " NUM NEW SWITCHES: " << numSwitchesUsed << " NUM HARDEN: " << numHardenUsed << " NUM FACILITY: " << numFacilityUsed <<  " OBJECTIVE: " << incumbentObjective << endl;

		if (CurrentObjective[0].val < -1e-4) {

			env = SubSolver[NumScenarios].cplex.getEnv();	

			isInfeasible = true;

			++bendersIter;
			cout << "Iteration #" << bendersIter << endl;

			int ScenarioIndex = CurrentObjective[0].ind;
			int numDamagedEdges;
			ScenarioIncluded[ScenarioIndex] = true;
			addedScenarios.push_back(ScenarioIndex);

			cout << "SubProblem #" << ScenarioIndex << " selected!" << endl;

			clonelineExistsConstraint.push_back(IloRangeArray(env));
			clonelineConstraint1.push_back(IloRangeArray(env));
			clonelineConstraint2.push_back(IloRangeArray(env));
			clonelineReactiveConstraint1.push_back(IloRangeArray(env));
			clonelineReactiveConstraint2.push_back(IloRangeArray(env));
			clonedirectionConstraint.push_back(IloRangeArray(env));
			cloneswitchConstraint1.push_back(IloRangeArray(env));
			cloneswitchConstraint2.push_back(IloRangeArray(env));
			cloneswitchReactiveConstraint1.push_back(IloRangeArray(env));
			cloneswitchReactiveConstraint2.push_back(IloRangeArray(env));
			clonevariationConstraint1.push_back(IloRangeArray(env));
			clonevariationConstraint2.push_back(IloRangeArray(env));
			clonevariationReactiveConstraint1.push_back(IloRangeArray(env));
			clonevariationReactiveConstraint2.push_back(IloRangeArray(env));
			clonedamageConstraint.push_back(IloRangeArray(env));
			cloneloadConstraint.push_back(IloRangeArray(env));
			cloneloadReactiveConstraint.push_back(IloRangeArray(env));
			clonegenerationConstraint.push_back(IloRangeArray(env));
			clonegenerationReactiveConstraint.push_back(IloRangeArray(env));
			clonebalanceConstraint.push_back(IloRangeArray(env));
			clonebalanceReactiveConstraint.push_back(IloRangeArray(env));
			clonemicrogridConstraint.push_back(IloRangeArray(env));
			clonecycleConstraint.push_back(IloRangeArray(env));
			clonecycleSubConstraint.push_back(IloRangeArray(env));
			cloneswitchableConstraint.push_back(IloRangeArray(env));
			clonecriticalServeConstraint.push_back(IloRangeArray(env));
			clonecriticalServeReactiveConstraint.push_back(IloRangeArray(env));
			cloneloadServeConstraint.push_back(IloRangeArray(env));
			cloneloadServeReactiveConstraint.push_back(IloRangeArray(env));
			cloneswitchImplicationConstraint1.push_back(IloRangeArray(env));
			//cloneswitchImplicationConstraint2.push_back(IloRangeArray(env));
			clonehardenConstraint.push_back(IloRangeArray(env));
			clonevoltageConstraint.push_back(IloRangeArray(env));
			clonevoltageOffsetConstraint.push_back(IloRangeArray(env));
			clonevoltageVariationConstraint1.push_back(IloRangeArray(env));
			clonevoltageVariationConstraint2.push_back(IloRangeArray(env));

			// VOLTAGE INFEASIBILITY CONSTRAINTS
			// clonevoltageInfeasibilityConstraint.push_back(IloRangeArray(env));

			clonelineUseVariable.push_back(IloNumVarArray(env));
			clonelineExistsVariable.push_back(IloNumVarArray(env));
			clonelineDirectionVariable.push_back(IloNumVarArray(env));
			clonelineCycleVariable.push_back(IloNumVarArray(env));
			clonelineHardenVariable.push_back(IloNumVarArray(env));
			cloneswitchUseVariable.push_back(IloNumVarArray(env));
			cloneswitchCycleVariable.push_back(IloNumVarArray(env));
			cloneflowVariable.push_back(IloNumVarArray(env));
			cloneflowReactiveVariable.push_back(IloNumVarArray(env));
			clonevoltageVariable.push_back(IloNumVarArray(env));
			clonevoltageOffsetVariable.push_back(IloNumVarArray(env));
			cloneloadVariable.push_back(IloNumVarArray(env));
			cloneloadReactiveVariable.push_back(IloNumVarArray(env));
			cloneloadServeVariable.push_back(IloNumVarArray(env));
			clonegeneratorVariable.push_back(IloNumVarArray(env));
			clonegeneratorReactiveVariable.push_back(IloNumVarArray(env));
			//clonemicrogridVariable.push_back(IloNumVarArray(env));
			clonefacilityVariable.push_back(IloNumVarArray(env));

			
			vector<double> redObjectiveVector(ColSize,0.0);

			PopulateScenarioData(
				SubSolver[NumScenarios],
				hashTableVertex,
				hashTableEdge,
				hashTableGenerator,
				hashTableLoad,
				NODES,
				EDGES,
				GENERATORS,
				LOADS,
				LINECODES,
				LINE_CONSTRUCTION_COST,
				MICROGRID_COST,
				MICROGRID_FIXED_COST,
				IS_CRITICAL_LOAD,
				MAX_MICROGRID,
				HARDEN_COST,
				LINE_SWITCH_COST,
				LoadMet,
				CriticalLoadMet,
				PhaseVariation,
				G,
				CYCLES,
				cloneSubObjective,
				redObjectiveVector,
				clonelineExistsConstraint[bendersIter], // CONSTRAINTS
				clonelineConstraint1[bendersIter], // CONSTRAINTS
				clonelineConstraint2[bendersIter],
				clonelineReactiveConstraint1[bendersIter], // CONSTRAINTS
				clonelineReactiveConstraint2[bendersIter],
				clonedirectionConstraint[bendersIter],
				cloneswitchConstraint1[bendersIter],
				cloneswitchConstraint2[bendersIter],
				cloneswitchReactiveConstraint1[bendersIter],
				cloneswitchReactiveConstraint2[bendersIter],
				clonevariationConstraint1[bendersIter],
				clonevariationConstraint2[bendersIter],
				clonevariationReactiveConstraint1[bendersIter],
				clonevariationReactiveConstraint2[bendersIter],
				clonedamageConstraint[bendersIter],
				cloneloadConstraint[bendersIter],
				cloneloadReactiveConstraint[bendersIter],
				clonegenerationConstraint[bendersIter],
				clonegenerationReactiveConstraint[bendersIter],
				clonebalanceConstraint[bendersIter],
				clonebalanceReactiveConstraint[bendersIter],
				clonemicrogridConstraint[bendersIter],
				clonecycleConstraint[bendersIter],
				clonecycleSubConstraint[bendersIter],
				cloneswitchableConstraint[bendersIter],
				clonecriticalServeConstraint[bendersIter],
				clonecriticalServeReactiveConstraint[bendersIter],
				cloneloadServeConstraint[bendersIter],
				cloneloadServeReactiveConstraint[bendersIter],
				cloneswitchImplicationConstraint1[bendersIter],
				//cloneswitchImplicationConstraint2[bendersIter],
				clonehardenConstraint[bendersIter],
				clonevoltageConstraint[bendersIter],
				clonevoltageOffsetConstraint[bendersIter],
				clonevoltageVariationConstraint1[bendersIter],
				clonevoltageVariationConstraint2[bendersIter],
				clonelineUseVariable[bendersIter], // VARIABLES
				clonelineExistsVariable[bendersIter], // VARIABLES
				clonelineDirectionVariable[bendersIter],
				clonelineCycleVariable[bendersIter],
				clonelineHardenVariable[bendersIter],
				cloneswitchUseVariable[bendersIter],
				cloneswitchCycleVariable[bendersIter],
				cloneflowVariable[bendersIter],
				cloneflowReactiveVariable[bendersIter],
				clonevoltageVariable[bendersIter],
				clonevoltageOffsetVariable[bendersIter],
				cloneloadVariable[bendersIter],
				cloneloadReactiveVariable[bendersIter],
				cloneloadServeVariable[bendersIter],
				clonegeneratorVariable[bendersIter],
				clonegeneratorReactiveVariable[bendersIter],
				//clonemicrogridVariable[bendersIter],
				clonefacilityVariable[bendersIter],
				numDamagedEdges,
				NumPhases,
				ScenarioIndex,
				NumScenarios,
				discreteMicrogrid,
				damage_str,
				root_str,
				HardenedDisabled,
				DISABLED[ScenarioIndex],
				HARDENED_DISABLED[ScenarioIndex],
				newImpedanceMultiplier,
				LDFIndicator);


		
			SubProblem[NumScenarios].add(clonelineExistsConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonelineConstraint1[bendersIter]);
			SubProblem[NumScenarios].add(clonelineConstraint2[bendersIter]);
			SubProblem[NumScenarios].add(clonelineReactiveConstraint1[bendersIter]);
			SubProblem[NumScenarios].add(clonelineReactiveConstraint2[bendersIter]);
			SubProblem[NumScenarios].add(clonedirectionConstraint[bendersIter]);
			SubProblem[NumScenarios].add(cloneswitchConstraint1[bendersIter]);
			SubProblem[NumScenarios].add(cloneswitchConstraint2[bendersIter]);
			SubProblem[NumScenarios].add(cloneswitchReactiveConstraint1[bendersIter]);
			SubProblem[NumScenarios].add(cloneswitchReactiveConstraint2[bendersIter]);
			SubProblem[NumScenarios].add(clonevariationConstraint1[bendersIter]);
			SubProblem[NumScenarios].add(clonevariationConstraint2[bendersIter]);
			SubProblem[NumScenarios].add(clonevariationReactiveConstraint1[bendersIter]);
			SubProblem[NumScenarios].add(clonevariationReactiveConstraint2[bendersIter]);
			SubProblem[NumScenarios].add(clonedamageConstraint[bendersIter]);
			SubProblem[NumScenarios].add(cloneloadConstraint[bendersIter]);
			SubProblem[NumScenarios].add(cloneloadReactiveConstraint[bendersIter]);
			//SubProblem[NumScenarios].add(clonegenerationConstraint[bendersIter]);
			//SubProblem[NumScenarios].add(clonegenerationReactiveConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonebalanceConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonebalanceReactiveConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonemicrogridConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonecycleConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonecycleSubConstraint[bendersIter]);
			SubProblem[NumScenarios].add(cloneswitchableConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonecriticalServeConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonecriticalServeReactiveConstraint[bendersIter]);
			SubProblem[NumScenarios].add(cloneloadServeConstraint[bendersIter]);
			SubProblem[NumScenarios].add(cloneloadServeReactiveConstraint[bendersIter]);
			SubProblem[NumScenarios].add(cloneswitchImplicationConstraint1[bendersIter]);
			//SubProblem[NumScenarios].add(cloneswitchImplicationConstraint2[bendersIter]);
			SubProblem[NumScenarios].add(clonehardenConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonevoltageConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonevoltageOffsetConstraint[bendersIter]);
			SubProblem[NumScenarios].add(clonevoltageVariationConstraint1[bendersIter]);
			SubProblem[NumScenarios].add(clonevoltageVariationConstraint2[bendersIter]);

			// VOLTAGE INFEASIBILITY CONSTRAINTS
			// SubProblem[NumScenarios].add(clonevoltageInfeasibilityConstraint[bendersIter]);

			SubProblem[NumScenarios].add(clonelineUseVariable[bendersIter]);
			SubProblem[NumScenarios].add(clonelineExistsVariable[bendersIter]);
			SubProblem[NumScenarios].add(clonelineDirectionVariable[bendersIter]);
			SubProblem[NumScenarios].add(clonelineCycleVariable[bendersIter]);
			SubProblem[NumScenarios].add(clonelineHardenVariable[bendersIter]);
			SubProblem[NumScenarios].add(cloneswitchUseVariable[bendersIter]);
			SubProblem[NumScenarios].add(cloneswitchCycleVariable[bendersIter]);
			SubProblem[NumScenarios].add(cloneflowVariable[bendersIter]);
			SubProblem[NumScenarios].add(cloneflowReactiveVariable[bendersIter]);
			SubProblem[NumScenarios].add(clonevoltageVariable[bendersIter]);
			SubProblem[NumScenarios].add(clonevoltageOffsetVariable[bendersIter]);
			SubProblem[NumScenarios].add(cloneloadVariable[bendersIter]);
			SubProblem[NumScenarios].add(cloneloadReactiveVariable[bendersIter]);
			SubProblem[NumScenarios].add(cloneloadServeVariable[bendersIter]);
			SubProblem[NumScenarios].add(clonegeneratorVariable[bendersIter]);
			SubProblem[NumScenarios].add(clonegeneratorReactiveVariable[bendersIter]);
			//SubProblem[NumScenarios].add(clonemicrogridVariable[bendersIter]);
			//SubProblem[NumScenarios].add(clonefacilityVariable[bendersIter]);


			cout << "SubProblem #" << ScenarioIndex << " added!" << endl;

			cout << "New constraint size: " << BaseLineVarConstraint.getSize() << endl;

			//cout << "LINES!" << endl;
			for (int j = 0; j < NumEdges; ++j) {
				if (DISABLED[ScenarioIndex][j].data == true) {
					BaseLineVarConstraint.add(IloRange(env, 0.0, IloInfinity));
				} else {
					BaseLineVarConstraint.add(IloRange(env, 0.0, 0.0));
				}
				BaseLineVarConstraint[bendersIter*NumEdges + j].setLinearCoef(lineUseVariable[NumScenarios][j], 1.0);
				BaseLineVarConstraint[bendersIter*NumEdges + j].setLinearCoef(clonelineUseVariable[bendersIter][j], -1.0);
				SubProblem[NumScenarios].add(BaseLineVarConstraint[bendersIter*NumEdges + j]);
			}
			//cout << "SWITCHES!" << endl;
			for (int j = 0; j < NumEdges; ++j) {
				BaseSwitchVarConstraint.add(IloRange(env, 0.0, IloInfinity));
				BaseSwitchVarConstraint[bendersIter*NumEdges + j].setLinearCoef(switchUseVariable[NumScenarios][j], 1.0);
				BaseSwitchVarConstraint[bendersIter*NumEdges + j].setLinearCoef(cloneswitchUseVariable[bendersIter][j], -1.0);
				SubProblem[NumScenarios].add(BaseSwitchVarConstraint[bendersIter*NumEdges + j]);
			}
			//cout << "HARDENS!" << endl;
			for (int j = 0; j < NumEdges; ++j) {
				if (HARDENED_DISABLED[ScenarioIndex][j].data == true) {
					BaseHardenVarConstraint.add(IloRange(env, 0.0, IloInfinity));
				} else {
					BaseHardenVarConstraint.add(IloRange(env, 0.0, 0.0));
				}
				BaseHardenVarConstraint[bendersIter*NumEdges + j].setLinearCoef(lineHardenVariable[NumScenarios][j], 1.0);
				BaseHardenVarConstraint[bendersIter*NumEdges + j].setLinearCoef(clonelineHardenVariable[bendersIter][j], -1.0);
				SubProblem[NumScenarios].add(BaseHardenVarConstraint[bendersIter*NumEdges + j]);
			}
			// ------------------------------------------------------------------------------------------
			// GENERATION CONSTRAINTS
			numCon = -1;
			int numPrevCon = BaseFacilityVarConstraint.getSize();
			for (int j = 0; j < NODES.size(); ++j) {
				for (int k = 0; k < NumPhases; ++k) {
					++numCon;
					BaseFacilityVarConstraint.add(IloRange(env, -MAXREALPHASE[j][k], IloInfinity));
					BaseFacilityVarReactiveConstraint.add(IloRange(env, -MAXREACTIVEPHASE[j][k], IloInfinity));
					BaseFacilityVarConstraint[numPrevCon + numCon].setLinearCoef(clonegeneratorVariable[bendersIter][j*NumPhases + k], -1.0);
					BaseFacilityVarReactiveConstraint[numPrevCon + numCon].setLinearCoef(clonegeneratorReactiveVariable[bendersIter][j*NumPhases + k], -1.0);
					for (int l = 0; l < GENERATORS.size(); ++l) {	
						if (hashTableVertex[GENERATORS[l].node_id] == j && GENERATORS[l].hasphase[k]) {
							BaseFacilityVarConstraint[numPrevCon + numCon].setLinearCoef(facilityVariable[NumScenarios][l], MAX_ADDED[l]);
							BaseFacilityVarReactiveConstraint[numPrevCon + numCon].setLinearCoef(facilityVariable[NumScenarios][l], MAX_ADDED[l]);
						}
					}
					SubProblem[NumScenarios].add(BaseFacilityVarConstraint[numPrevCon + numCon]);	
					SubProblem[NumScenarios].add(BaseFacilityVarReactiveConstraint[numPrevCon + numCon]);	
				}
			}

			cout << "Finished modifications!" << endl; 
		} else {
			if (FirstVoltageCut) {
				CPUTIME_SBD = (double) (clock() - begin) / (double) CLOCKS_PER_SEC;
				FirstVoltageCut = false;
			}
			// TERMINATE IF MAX_CPUTIME REACHED
			if (((double) (clock() - begin) / (double) CLOCKS_PER_SEC) > MAX_CPUTIME-60.0) {
				break;
			}


			cout << "CPUTIME: " << (double) (clock() - begin) / (double) CLOCKS_PER_SEC << endl;
			vector<vector<vector<double> > > voltageFeasible(
				NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));
			for (int i = 0; i < NumScenarios; ++i) {
				for (int j = 0; j < NumNodes; ++j) {
					for (int k = 0; k < NumPhases; ++k) {
						voltageFeasible[i][j][k] = 24900; //NODES[j].refVoltage;
					}
				}
				voltageInfeasibleScenarios[i] = false;
			}


			/*cloneFlowVals.clear();
			cloneFlowVals.push_back(IloNumArray(env));
			SubSolver[NumScenarios].cplex.getValues(cloneFlowVals[0], cloneflowVariable[0]);


			cout << "MASTER subxf: " << lineusevals[hashTableEdge["subxf"]] << " " << lineswitchvals[hashTableEdge["subxf"]] << endl;
			cout << "LINE subxf: " << cloneLineUseVals[0][hashTableEdge["subxf"]] << " " << cloneSwitchUseVals[0][hashTableEdge["subxf"]] << endl;


			cout << cloneFlowVals[0][hashTableEdge["subxf"]*NumPhases] << " " 
				<< cloneFlowVals[0][hashTableEdge["subxf"]*NumPhases + 1] << " " 
				<< cloneFlowVals[0][hashTableEdge["subxf"]*NumPhases + 2] << endl;*/

			/*for (int j = 0; j < NumEdges; ++j) {
				cout << cloneLineExistsVals[0][j] << " ";
			}
			cout << endl;

			for (int j = 0; j < NumEdges; ++j) {
				cout << (cloneLineUseVals[0][j] - cloneSwitchUseVals[0][j]) << " ";
			}
			cout << endl;*/

			vector<vector<bool> > activeGenerators(NumScenarios, vector<bool>(NumNodes, false));
			vector<vector<bool> > activeLoads(NumScenarios, vector<bool>(NumNodes, false));
			int numPrevComp = scenarioList.size();
			string DN_type = extent_str;
			DN_type[0] = tolower(DN_type[0]);



			vector<vector<int> > allForests(NumScenarios, vector<int>(NumEdges, 0));
			vector<vector<int> > allLoads(NumScenarios, vector<int>(NumLoads, 0));
			vector<vector<vector<double> > > allVoltage(NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));
			vector<vector<vector<double> > > allRealGen(NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));
			vector<vector<vector<double> > > allReactiveGen(NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));
			vector<vector<vector<double> > > allRealLoad(NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));
			vector<vector<vector<double> > > allReactiveLoad(NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));

			string gridlabd_path;
			if (ExactSolver) {
				gridlabd_path = "gridlabd_output/Exact/" + string(LDFIndicator ? "LDF/" : "MCF/") + root_str + damage_str;
			} else {
				gridlabd_path = "gridlabd_output/BVNDS/" + string(LDFIndicator ? "LDF/" : "MCF/") + root_str + damage_str;
			}
			cout << gridlabd_path << endl;
			populateGridLabDInput(
					NODES,
					EDGES,
					GENERATORS,
					LOADS,
					LINECODES,
					hashTableVertex,
					hashTableEdge,
					hashTableGenerator,
					hashTableLoad,
					cloneLineUseVals,
					cloneLineDirectionVals,
					cloneSwitchUseVals,
					cloneGeneratorVals,
					cloneGeneratorReactiveVals,
					cloneLoadVals,
					cloneLoadReactiveVals,
					cloneLoadServeVals,
					cloneVoltageVals,
					subLineUseVals,
					subLineDirectionVals,
					subSwitchUseVals,
					subGeneratorVals,
					subGeneratorReactiveVals,
					subLoadVals,
					subLoadReactiveVals,
					subLoadServeVals,
					subVoltageVals,
					ScenarioIncluded,
					addedScenarios,
					NumScenarios,
					NumNodes,
					NumEdges,
					NumLoads,
					NumPhases,
					allForests,
					allLoads,
					allVoltage,
					allRealGen,
					allReactiveGen,
					allRealLoad,
					allReactiveLoad,
					activeGenerators,
					activeLoads);



			for (int i = 0; i < NumScenarios; ++i) {
				vector<int> newComponent(NumEdges + NumLoads + NumGenerators, 0);
				for (int j = 0; j < NumEdges; ++j) { 
					/*int edgeVal = 0;
					if (allForests[i][j] != 0) {
						edgeVal = 1;
					}*/
					newComponent[j] = allForests[i][j];
				}
				for (int j = 0; j < NumLoads; ++j) {
					/*int loadVal = 0;
					if (activeLoads[i][j]) {
						loadVal = 1;
					}*/
					newComponent[NumEdges+j] = allLoads[i][j];
				}
				for (int j = 0; j < NumGenerators; ++j) {
					int generatorVal = 0;
					if (fabs((double)facilityvals[j] - 1.0) < 1e-4) {
						generatorVal = 1;
					}
					newComponent[NumEdges+NumLoads+j] = generatorVal;
				}
				set<vector<int> >::iterator it = componentSolutions[i].find(newComponent);
				if (it == componentSolutions[i].end()) {
					componentSolutions[i].insert(newComponent);
				} else {
					cout << "********* PREVIOUS COMPONENT FOUND!!!! **********" << endl;
				}
			}		

			cout << "OUTPUTTING TO GRIDLABDPATH" << endl;
			// SORT ORDERS
			for (int ComponentOrder = 1; ComponentOrder < 2; ++ComponentOrder) {
				GridLabDInput(
					NODES,
					EDGES,
					GENERATORS,
					LOADS,
					LINECODES,
					hashTableVertex,
					hashTableEdge,
					hashTableGenerator,
					hashTableLoad,
					ScenarioIncluded,
					addedScenarios,
					NumScenarios,
					NumNodes,
					NumEdges,
					NumPhases,
					allForests,
					allVoltage,
					allRealGen,
					allReactiveGen,
					allRealLoad,
					allReactiveLoad,
					ComponentOrder,
                                        ComponentOrder,
					gridlabd_path);

				// PREVIOUS GRIDLABD CALL
				/*if (checkGridLabD(NumScenarios, newImpedanceMultiplier, 
						hashTableVertex, voltageFeasible))
					cout << "FEASIBLE" << endl;
				else 
					cout << "INFEASIBLE" << endl;*/

				/*checkGridLabD_v2(NumScenarios, newImpedanceMultiplier,
					gridlabd_path, DN_type, ComponentOrder);


				readGridLabDOutput(
					NODES,
					EDGES,
					GENERATORS,
					LOADS,
					LINECODES,
					hashTableVertex,
					hashTableEdge,
					hashTableGenerator,
					hashTableLoad,	
					scenarioList,
					nodesList,
					infNodeList,
					edgesList,
					NumScenarios,
					NumNodes,
					NumEdges,
					NumPhases,
					allForests,
					activeGenerators,
					gridlabd_path);*/
			}
			break;

			//if (numVoltageCuts == 1) return;

			/*for (int l = numPrevComp; l < scenarioList.size(); ++l) {
				int scenarioInd = scenarioList[l];
				voltageInfeasibleScenarios[scenarioInd] = true;
			}


			cout << "NUMBER OF NEW COMPONENTS: " << (scenarioList.size() - numPrevComp) << endl;
			cout << "Number of voltage cuts: " << numVoltageCuts << " MAX " << MaxVoltageRounds << endl;
			if ((scenarioList.size() - numPrevComp) > 0 && numVoltageCuts < MaxVoltageRounds) {

				// GENERATE MORE CUTS USING RANDOM ORDERING
				for (int l = 0; l < 0; ++l) {
					GridLabDInput(
						NODES,
						EDGES,
						GENERATORS,
						LOADS,
						LINECODES,
						hashTableVertex,
						hashTableEdge,
						hashTableGenerator,
						hashTableLoad,
						ScenarioIncluded,
						addedScenarios,
						NumScenarios,
						NumNodes,
						NumEdges,
						NumPhases,
						allForests,
						allVoltage,
						allRealGen,
						allReactiveGen,
						allRealLoad,
						allReactiveLoad,
						2,
						l+2,
						gridlabd_path);


					checkGridLabD_v2(NumScenarios, newImpedanceMultiplier,
						gridlabd_path, DN_type, l+2);


					readGridLabDOutput(
						NODES,
						EDGES,
						GENERATORS,
						LOADS,
						LINECODES,
						hashTableVertex,
						hashTableEdge,
						hashTableGenerator,
						hashTableLoad,	
						scenarioList,
						nodesList,
						infNodeList,
						edgesList,
						NumScenarios,
						NumNodes,
						NumEdges,
						NumPhases,
						allForests,
						activeGenerators,
						gridlabd_path);
				}




				//env = SubSolver[NumScenarios].cplex.getEnv();
				//!checkGridLabD(NumScenarios, newImpedanceMultiplier, 
				//hashTableVertex, voltageFeasible)) {

				cout << "VOLTAGE INFEASIBLE.. REINITIALIZING" << endl;
				numVoltageCuts++;			
	
				//SubSolver[NumScenarios].cplex.exportModel("LinDistFlow.lp");
				for (int i = 0; i < bendersIter+1; ++i) {
					//numInfeasibilityCuts = clonevoltageInfeasibilityConstraint[i].getSize() - 1;
					isInfeasible = true;

					//SubProblem[NumScenarios].remove(clonevoltageInfeasibilityConstraint[i]);
					//clonevoltageInfeasibilityConstraint[i].end();
					//clonevoltageInfeasibilityConstraint[i] = IloRangeArray(env);
					//SubProblem[NumScenarios].add(clonevoltageInfeasibilityConstraint[i]);

					for (int l = numPrevComp; l < scenarioList.size(); ++l) {
						int scenarioInd = scenarioList[l];

						double CurrRHS = 1.0;
						numInfeasibilityCuts++;
						InfeasibilityCuts.add(IloRange(env, CurrRHS, IloInfinity));


						vector<bool> activeNodes(NumNodes, false);
						vector<bool> activeEdges(NumEdges, false);
						//vector<bool> activeGenerators(NumNodes, false);
						for (list<string>::iterator it = edgesList[l].begin(); it != edgesList[l].end(); ++it) {
							//cout << *it << " ";
							activeEdges[hashTableEdge[*it]] = true;
							//activeGenerators[hashTableVertex[EDGES[hashTableEdge[*it]].node1id]] = true;
							//activeGenerators[hashTableVertex[EDGES[hashTableEdge[*it]].node2id]] = true;
							//activeNodes[hashTableVertex[EDGES[hashTableEdge[*it]].node1id]] = true;
							//activeNodes[hashTableVertex[EDGES[hashTableEdge[*it]].node2id]] = true;
						}			
						//cout << endl;		
						for (list<string>::iterator it = nodesList[l].begin(); it != nodesList[l].end(); ++it) {
							activeNodes[hashTableVertex[*it]] = true;
						}

						//cout << "ADDING CUT #" << l << endl;

						// cout << "NumInfCuts: " << numInfeasibilityCuts << endl;
						// L-SHAPED TOPOLOGY CUT
						for (int j = 0; j < NumEdges; ++j) {
							if (activeEdges[j]) {// && (EDGES[j].id[0] == 'n' || EDGES[j].id[0] == 'o')) {
								if (allForests[scenarioInd][j] == 0) {
									cout << "(" << EDGES[j].id << "," << allForests[i][j] << "," << cloneLineExistsVals[i][j] << ") " << endl;
									//InfeasibilityCuts[numInfeasibilityCuts].setLinearCoef(clonelineExistsVariable[i][j],1.0);
								} else {
									//cout << EDGES[j].id << " ";
									InfeasibilityCuts[numInfeasibilityCuts].setLinearCoef(clonelineExistsVariable[i][j],-1.0);
									CurrRHS -= 1.0;
								}
							}
						}
						for (int j = 0; j < LOADS.size(); ++j) {
							if (activeNodes[hashTableVertex[LOADS[j].node_id]] && infNodeList[l] == LOADS[j].node_id) {
								//cout << "LOAD " << LOADS[j].id << " NODE " << LOADS[j].node_id << endl;
								if (allLoads[scenarioInd][j] == 0) {
									//cout << LOADS[j].id << " ";
									//InfeasibilityCuts[numInfeasibilityCuts].setLinearCoef(cloneloadServeVariable[i][j],1.0);
								} else {
									//cout << "-" << LOADS[j].id << " ";
									InfeasibilityCuts[numInfeasibilityCuts].setLinearCoef(cloneloadServeVariable[i][j],-1.0);
									CurrRHS -= 1.0;
								}
							}
						}
						for (int j = 0; j < NumGenerators; ++j) {
							if (GENERATORS[j].id[0] != 's' && activeNodes[hashTableVertex[GENERATORS[j].node_id]]) {
								if (fabs((double)facilityvals[j] - 0.0) < 1e-4) { //   !activeGenerators[scenarioInd][hashTableVertex[GENERATORS[j].node_id]]) {
									InfeasibilityCuts[numInfeasibilityCuts].setLinearCoef(facilityVariable[NumScenarios][j],1.0);
								} else {
									//InfeasibilityCuts[numInfeasibilityCuts].setLinearCoef(facilityVariable[NumScenarios][j],-1.0);
									//CurrRHS -= 1.0;
								}
							}
						}
						//cout << endl;

						InfeasibilityCuts[numInfeasibilityCuts].setLB(CurrRHS);
						SubProblem[NumScenarios].add(InfeasibilityCuts[numInfeasibilityCuts]);
					}

				}
			
			} else {
				cout << "VOLTAGE FEASIBLE SOLUTION" << endl;
			}*/
		}


	} while (isInfeasible);
	time_t end = clock();

	cout << "TIMED OUT: " << (TIMED_WITHOUT_SOLUTION ? "LIMIT REACHED" : "SOLVED") << endl;
	cout << "STATUS: Feasible OBJECTIVE: " << incumbentObjective << endl;
	cout << "MINIMUM FEASIBLE OBJECTIVE: " << MIN_FEAS_OBJ << endl;
	cout << "NUM SCENARIOS GENERATED: " << bendersIter+1 << endl;
	cout << "CPUTIME_SBD: " << CPUTIME_SBD << endl;
	cout << "CPUTIME_GLD: " << (double) (end - begin) / (double) CLOCKS_PER_SEC << endl;

	cout << "NUMBER OF VOLTAGE ROUNDS: " << numVoltageCuts << endl;
	cout << "NUMBER OF VOLTAGE CUTS: " << numInfeasibilityCuts << endl;

	int infeasibleScenarioSum = 0;
	for (int i = 0; i < NumScenarios; ++i) {
		if (voltageInfeasibleScenarios[i]) {
			infeasibleScenarioSum++;
		}
	}
	cout << "NUMBER OF VOLTAGE INFEASIBLE SCENARIOS: " << infeasibleScenarioSum << endl;

	if (CHECK_RESULTS) {
		env = SubSolver[NumScenarios].cplex.getEnv();
		IloNumArray lineusevals(env);
		IloNumArray lineswitchvals(env);
		IloNumArray linehardenvals(env);
		//IloNumArray microgridvals(env);				
		IloNumArray facilityvals(env);

		for (int j = 0; j < NumEdges; ++j) {
			lineusevals.add(incumbentSolution[j]);
			lineswitchvals.add(incumbentSolution[j + NumEdges]);
			linehardenvals.add(incumbentSolution[j + 2*NumEdges]);
		}
		for (int j = 0; j < NumGenerators; ++j) {
			/*for (int k = 0; k < NumPhases; ++k) {
				microgridvals.add(incumbentSolution[j*NumPhases + k + 3*NumEdges]);
			}*/
			facilityvals.add(incumbentSolution[j + 3*NumEdges]);
		}

		std::stringstream ss1;
		ss1 << LoadMet;
		std::stringstream ss2;
		ss2 << CriticalLoadMet;
		std::stringstream ss3;
		ss3 << discreteMicrogrid;
		string critical_str = ss1.str() + "_" + ss2.str() + "_" + ss3.str();

		string str;
		if (ExactSolver) {
			str = "results/Exact/" + string(LDFIndicator ? "LDF/" : "MCF/") + root_str + "solution_" + critical_str + "_" + damage_str + ".txt";
		} else {
			str = "results/BVNDS/" + string(LDFIndicator ? "LDF/" : "MCF/") + root_str + "solution_" + critical_str + "_" + damage_str + ".txt";
		}
		cout << str << endl;
		ofstream outputSolution(str.c_str());

		outputSolution << "NODES: " << endl;
		for (int j = 0; j < NODES.size(); ++j) {
			outputSolution << "NODE " << NODES[j].id << " (" << NODES[j].x << "," << NODES[j].y << ")" << endl;
		}

		outputSolution << "FLOW VALUES: " << endl;
		for (int j = 0; j < EDGES.size(); ++j) {
			outputSolution << "EDGE " << EDGES[j].id 
				 << " (" << (EDGES[j].hasphase[0] == true ? 1 : 0) << "," << (EDGES[j].hasphase[1] == true ? 1 : 0) << "," << (EDGES[j].hasphase[2] == true ? 1 : 0) << ") " << EDGES[j].capacity
				 << " (" << EDGES[j].node1id << "," << EDGES[j].node2id << ") (" << lineusevals[j] << "," << lineswitchvals[j] << "," << linehardenvals[j] << ")" << endl;
		}
		outputSolution << "GENERATION VALUES: " << endl;
		for (int j = 0; j < GENERATORS.size(); ++j) {
			outputSolution << "GENERATOR " << GENERATORS[j].id 
				 << " (" << (GENERATORS[j].hasphase[0] == true ? 1 : 0) << "," << (GENERATORS[j].hasphase[1] == true ? 1 : 0) << "," << (GENERATORS[j].hasphase[2] == true ? 1 : 0) << ") "
				 << MAX_ADDED[j] << " (" << facilityvals[j] << ")" << endl;
			
		}

		// FIND ALL FORESTS:
		vector<vector<int> > allForests(NumScenarios, vector<int>(NumEdges, 0));
		vector<vector<vector<double> > > allRealGen(NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));
		vector<vector<vector<double> > > allReactiveGen(NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));
		vector<vector<vector<double> > > allRealLoad(NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));
		vector<vector<vector<double> > > allReactiveLoad(NumScenarios, vector<vector<double> >(NumNodes, vector<double>(NumPhases, 0.0)));
		cout << "Number of included subproblems: " << cloneLineUseVals.size() << endl;
		for (int j = 0; j < cloneLineUseVals.size(); ++j) {
			for (int e = 0; e < NumEdges; ++e) {
				if (fabs((double)cloneLineUseVals[j][e] - 1.0) < 1e-4 && fabs((double)cloneSwitchUseVals[j][e] - 0.0) < 1e-4) {
					if (fabs((double) cloneLineDirectionVals[j][2*e] - 1.0) < 1e-4) {
						allForests[addedScenarios[j]][e] = -1;
					} else {
						allForests[addedScenarios[j]][e] = 1;

					}
					//cout << "1 ";
				} else {
					//cout << "0 ";
				}
			}
			for (int l = 0; l < NumNodes; ++l) {
				for (int k = 0; k < NumPhases; ++k) {
					allRealGen[addedScenarios[j]][l][k] = cloneGeneratorVals[j][l*NumPhases + k];
					allRealLoad[addedScenarios[j]][l][k] = cloneLoadVals[j][l*NumPhases + k];
					allReactiveLoad[addedScenarios[j]][l][k] = cloneLoadReactiveVals[j][l*NumPhases + k];
				}
			} 
		}
		for (int j = 0; j < NumScenarios; ++j) {
			if (!ScenarioIncluded[j]) {
				for (int e = 0; e < NumEdges; ++e) {
					if (fabs((double)subLineUseVals[j][e] - 1.0) < 1e-4 && fabs((double)subSwitchUseVals[j][e] - 0.0) < 1e-4) {
						if (fabs((double) subLineDirectionVals[j][2*e] - 1.0) < 1e-4) {
							allForests[j][e] = -1;
						} else {
							allForests[j][e] = 1;
						}
					}
				} 
				for (int l = 0; l < NumNodes; ++l) {
					for (int k = 0; k < NumPhases; ++k) {
						allRealGen[j][l][k] = subGeneratorVals[j][l*NumPhases + k];
						allRealLoad[j][l][k] = subLoadVals[j][l*NumPhases + k];
						allReactiveLoad[j][l][k] = subLoadReactiveVals[j][l*NumPhases + k];
					}
				}
			}
		}	
		//cout << "NUMBER OF DISTINCT TREES: " << allForests.size() << endl;
		/*ofstream out1("HarshaOutput/forests.txt");
		for (int e = 0; e < NumEdges; ++e) {
			out1 << EDGES[e].id << ",";
			for (int i = 0; i < allForests.size(); ++i) {
				out1 << allForests[i][e] << ",";
			}
			out1 << endl;
		}
		out1.close();*/
	

		outputSolution << "STATUS: " << STATUS << " OBJECTIVE: " << incumbentObjective << endl;
		outputSolution << "NUM SCENARIOS GENERATED: " << bendersIter+1 << endl;
		outputSolution << "CPUTIME_SBD: " << CPUTIME_SBD << endl;
		outputSolution << "CPUTIME_GLD: " << (double) (end - begin) / (double) CLOCKS_PER_SEC << endl;

		outputSolution << "NUMBER OF VOLTAGE ROUNDS: " << numVoltageCuts << endl;
		outputSolution << "NUMBER OF VOLTAGE CUTS: " << numInfeasibilityCuts << endl;
		outputSolution << "NUMBER OF VOLTAGE INFEASIBLE SCENARIOS: " << infeasibleScenarioSum << endl;

		// CHECK IF THERE ARE CYCLES IN THE SOLUTION:
		for (int i = 0; i < NumScenarios; ++i) {
			// CREATE GRAPH
			vector<graph_vertex> GSub(NODES.size());
			for (int j = 0; j < EDGES.size(); ++j) {
				if (allForests[i][j] != 0) {
					GSub[hashTableVertex[EDGES[j].node1id]].AdjList.push_back(hashTableVertex[EDGES[j].node2id]);
					GSub[hashTableVertex[EDGES[j].node1id]].EdgeID.push_back(EDGES[j].id);
					GSub[hashTableVertex[EDGES[j].node2id]].AdjList.push_back(hashTableVertex[EDGES[j].node1id]);
					GSub[hashTableVertex[EDGES[j].node2id]].EdgeID.push_back(EDGES[j].id);
				}
			}

			// CREATE CYCLE LIST AND FIND ALL CYCLES
			vector<vector<int> > CYCLESSub;
			detectCycles(GSub, CYCLESSub);
			cout << "SUBGRAPH " << i << " HAS " << CYCLESSub.size() << " CYCLES" << endl;
			list<vector<int> > components;
			findComponents(GSub, components);

			for (list<vector<int> >::iterator it = components.begin();
				it != components.end(); ++it) {
				int numGenerators = 0;

				vector<double> componentRealGen(NumPhases, 0.0);				
				vector<double> componentReactiveLoad(NumPhases, 0.0);
				cout << "NEW COMPONENT: SIZE = " << it->size() << ", ";
				for (int j = 0; j < it->size(); ++j) {
					cout << NODES[(*it)[j]].id << " ";
					for (int k = 0; k < NumPhases; ++k) {
						componentRealGen[k] += allRealGen[i][(*it)[j]][k];
						componentReactiveLoad[k] += allReactiveLoad[i][(*it)[j]][k];
					}		
				}
				cout << endl;
				for (int j = 0; j < it->size(); ++j) {
					for (int k = 0; k < NumPhases; ++k) {
						if (componentRealGen[k] > 1e-4) {
							allReactiveGen[i][(*it)[j]][k] = (allRealGen[i][(*it)[j]][k] / componentRealGen[k]) * componentReactiveLoad[k];
						}
					}
				}
			}
		}

		/*ofstream out2("HarshaOutput/loadsReal.txt");
		for (int e = 0; e < NumNodes; ++e) {
			out2 << NODES[e].id << ",";
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					out2 << allRealLoad[i][e][k] << ",";
				}
			}
			out2 << endl;
		}
		out2.close();
		ofstream out3("HarshaOutput/loadsReactive.txt");
		for (int e = 0; e < NumNodes; ++e) {
			out3 << NODES[e].id << ",";
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					out3 << allReactiveLoad[i][e][k] << ",";
				}
			}
			out3 << endl;
		}
		out3.close();
		ofstream out4("HarshaOutput/generatorsReal.txt");
		for (int e = 0; e < NumNodes; ++e) {
			out4 << NODES[e].id << ",";
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					out4 << allRealGen[i][e][k] << ",";
				}
			}
			out4 << endl;
		}
		out4.close();
		ofstream out5("HarshaOutput/generatorsReactive.txt");
		for (int e = 0; e < NumNodes; ++e) {
			out5 << NODES[e].id << ",";
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < NumPhases; ++k) {
					out5 << allReactiveGen[i][e][k] << ",";
				}
			}
			out5 << endl;
		}
		out5.close();
		ofstream out6("HarshaOutput/lines.txt");
		for (int e = 0; e < NumEdges; ++e) {
			out6 << EDGES[e].id << "," << EDGES[e].node1id << "," << EDGES[e].node2id
			     << "," << EDGES[e].length << "," << EDGES[e].capacity;
			for (int k = 0; k < NumPhases; ++k) {
				if (EDGES[e].hasphase[k]) {
					out6 << ",1";
				} else {
					out6 << ",0";
				}
			}
			if (!EDGES[e].istransformer) {
				map<int, lineCodeData>::iterator it = LINECODES.find(EDGES[e].linecode);

				out6 << "," << it->first;
				for (int k = 0; k < EDGES[e].NumPhases; ++k) {
					for (int l = 0; l < EDGES[e].NumPhases; ++l) {
						out6 << "," << (it->second.rmatrix[0][k][l]);
					}

				}
				for (int k = 0; k < EDGES[e].NumPhases; ++k) {
					for (int l = 0; l < EDGES[e].NumPhases; ++l) {
						out6 << "," << (it->second.xmatrix[0][k][l]);
					}

				}
			}
			out6 << endl;
		}
		out6.close();*/

		outputSolution << "MINIMUM FEASIBLE OBJECTIVE: " << MIN_FEAS_OBJ << endl;
		/*outputSolution << "BENDERS ITERATIONS:" << endl;
		for (int j = 0; j < IterationObjective.size(); ++j) {
			outputSolution << IterationObjective[j] << " " << IterationCPUTIME[j] << endl;
		}*/
	}

}

// FOLLOWING METHOD PERFORMS:
// 1. READ INITIAL PROBLEM DATA
// 2. SOLVE PROBLEM SEPARATELY
// 3. SOLVE MASTER PROBLEM USING BENDERS VARIABLE NEIGHBORHOOD DECOMPOSITION SEARCH
void generateInitialProblem (
	IloModel &MasterProblem,
	vector<IloModel> &SubProblem,
	Solver &MasterSolver,
	vector<Solver> &SubSolver,
	IloObjective &MasterObjective, IloRangeArray &MasterConstraints, 
	vector<IloObjective> &SubObjective, // OBJECTIVE 
	vector<IloRangeArray> &lineExistsConstraint, // CONSTRAINTS
	vector<IloRangeArray> &lineConstraint1, // CONSTRAINTS
	vector<IloRangeArray> &lineConstraint2,
	vector<IloRangeArray> &lineReactiveConstraint1,
	vector<IloRangeArray> &lineReactiveConstraint2,
	vector<IloRangeArray> &directionConstraint,
	vector<IloRangeArray> &switchConstraint1,
	vector<IloRangeArray> &switchConstraint2,
	vector<IloRangeArray> &switchReactiveConstraint1,
	vector<IloRangeArray> &switchReactiveConstraint2,
	vector<IloRangeArray> &variationConstraint1,
	vector<IloRangeArray> &variationConstraint2,
	vector<IloRangeArray> &variationReactiveConstraint1,
	vector<IloRangeArray> &variationReactiveConstraint2,
	vector<IloRangeArray> &damageConstraint,
	vector<IloRangeArray> &loadConstraint,
	vector<IloRangeArray> &loadReactiveConstraint,
	vector<IloRangeArray> &generationConstraint,
	vector<IloRangeArray> &generationReactiveConstraint,
	vector<IloRangeArray> &balanceConstraint,
	vector<IloRangeArray> &balanceReactiveConstraint,
	vector<IloRangeArray> &microgridConstraint,
	vector<IloRangeArray> &cycleConstraint,
	vector<IloRangeArray> &cycleSubConstraint,
	vector<IloRangeArray> &switchableConstraint,
	vector<IloRangeArray> &criticalServeConstraint,
	vector<IloRangeArray> &criticalServeReactiveConstraint,
	vector<IloRangeArray> &loadServeConstraint,
	vector<IloRangeArray> &loadServeReactiveConstraint,
	vector<IloRangeArray> &switchImplicationConstraint1,
	//vector<IloRangeArray> &switchImplicationConstraint2,
	vector<IloRangeArray> &hardenConstraint,
	vector<IloRangeArray> &voltageConstraint,
	vector<IloRangeArray> &voltageOffsetConstraint,
	vector<IloRangeArray> &voltageVariationConstraint1,
	vector<IloRangeArray> &voltageVariationConstraint2,
	vector<IloNumVarArray> &lineUseVariable, // VARIABLES
	vector<IloNumVarArray> &lineExistsVariable, // VARIABLES
	vector<IloNumVarArray> &lineDirectionVariable,
	vector<IloNumVarArray> &lineCycleVariable,
	vector<IloNumVarArray> &lineHardenVariable,
	vector<IloNumVarArray> &switchUseVariable,
	vector<IloNumVarArray> &switchCycleVariable,
	vector<IloNumVarArray> &flowVariable,
	vector<IloNumVarArray> &flowReactiveVariable,
	vector<IloNumVarArray> &voltageVariable,
	vector<IloNumVarArray> &voltageOffsetVariable,
	vector<IloNumVarArray> &loadVariable,
	vector<IloNumVarArray> &loadReactiveVariable,
	vector<IloNumVarArray> &loadServeVariable,
	vector<IloNumVarArray> &generatorVariable,
	vector<IloNumVarArray> &generatorReactiveVariable,
	//vector<IloNumVarArray> &microgridVariable,
	vector<IloNumVarArray> &facilityVariable,					
	int NumScenarios,
	double CriticalLoadMet,
	double discreteMicrogrid,
	string damage_str,
	string extent_str,
	string root_str,
	bool HardenedDisabled,
	bool ExactSolver,
	bool LDFIndicator,
	int MaxVoltageRounds,
	bool CHECK_RESULTS,
	double newImpedanceMultiplier) {

	bool SOLVE = true;
	bool addMoreGenerators = false;

	int NumPhases = 3;

	IloEnv env = MasterSolver.cplex.getEnv();
	//MasterSolver.cplex.extract(MasterProblem);
	MasterSolver.cplex.setOut(env.getNullStream());
	for (int i = 0; i < NumScenarios + 1; ++i) {
		IloEnv env = SubSolver[i].cplex.getEnv();
		//SubSolver[i].cplex.extract(SubProblem[i]);
		SubSolver[i].cplex.setOut(env.getNullStream());
	}

	vector<vector<vector<double> > > SubColumn (NumScenarios + 1, vector<vector<double> >(0));

	// GRAPH DATA
	vector<nodeData> NODES;
	vector<edgeData> EDGES;
	vector<generatorData> GENERATORS_mod;
	vector<loadData> LOADS;
	map<int, lineCodeData> LINECODES;

	// Large Expansion/Large Extent/balancedData/large-extent.txt for ALL transformers balanced (even separate lines)
	string graph_file = root_str + extent_str + ".txt";
	string transformer_file = root_str + "Transformer.DSS";
	string line_file = root_str + "Line.DSS";
	string linecode_file = root_str + "LineCode.DSS";
	cout << graph_file << endl;
	cout << transformer_file << endl;
	// string root_str = "Large Expansion/Large Extent/"

	// READ GRAPH DATA
	readTextData(graph_file, transformer_file, line_file, linecode_file, NODES, 
			GENERATORS_mod, LOADS, EDGES, LINECODES);

	// -- ADD GENERATOR TO EVERY NODE ------------------------------------------
	int numInitialGenerators = GENERATORS_mod.size();
	if (addMoreGenerators) {
		for (int i = 0; i < NODES.size(); ++i) {
			GENERATORS_mod.push_back(generatorData());
			GENERATORS_mod[numInitialGenerators + i].id = "gA" + NODES[i].id;
			GENERATORS_mod[numInitialGenerators + i].node_id = NODES[i].id;
			for (int k = 0; k < NumPhases; ++k) {
				GENERATORS_mod[numInitialGenerators + i].hasphase.push_back(true);
				GENERATORS_mod[numInitialGenerators + i].maxrealphase.push_back(0.0);
				GENERATORS_mod[numInitialGenerators + i].maxreactivephase.push_back(0.0);
			}
		}
	}
	// -------------------------------------------------------------------------

	// CREATE GRAPH AS DEFINED IN detect_cycles.hpp
	// Following Hash Table will be used to map node ids into integers for instant (O(1)) access to the graph.
	map<string, int> hashTableVertex;
	for (int i = 0; i < NODES.size(); ++i) {
		hashTableVertex.insert(pair<string, int>(NODES[i].id, i));
	}
	map<string, int> hashTableEdge;
	for (int i = 0; i < EDGES.size(); ++i) {
		hashTableEdge.insert(pair<string, int>(EDGES[i].id, i));
	}
	map<string, int> hashTableLoad;
	for (int i = 0; i < LOADS.size(); ++i) {
		hashTableLoad.insert(pair<string, int>(LOADS[i].id, i));
	}
	map<string, int> hashTableGenerator;
	for (int i = 0; i < GENERATORS_mod.size(); ++i) {
		cout << GENERATORS_mod[i].id << " " << i << endl;
		hashTableGenerator.insert(pair<string, int>(GENERATORS_mod[i].id, i));
	}

	// -- CALCULATE LINE LENGTHS ----------------------------------------------
	for (int i = 0; i < EDGES.size(); ++i) {
		double xdist = NODES[hashTableVertex[EDGES[i].node1id]].x - NODES[hashTableVertex[EDGES[i].node2id]].x;	
		double ydist = NODES[hashTableVertex[EDGES[i].node1id]].y - NODES[hashTableVertex[EDGES[i].node2id]].y;
		double edgeLength = pow(pow(xdist,2.0) + pow(ydist,2.0),0.5) * 69.0;
		EDGES[i].length = edgeLength;	
	}

	/*ofstream out("edges.txt");
	for (int i = 0; i < EDGES.size(); ++i) {
		out << EDGES[i].id << "," << EDGES[i].node1id << "," << EDGES[i].node2id << ","
		    << EDGES[i].length << endl;
	}
	out.close();*/

	// CREATE GRAPH
	vector<graph_vertex> G(NODES.size());
	for (int i = 0; i < EDGES.size(); ++i) {
		G[hashTableVertex[EDGES[i].node1id]].AdjList.push_back(hashTableVertex[EDGES[i].node2id]);
		G[hashTableVertex[EDGES[i].node1id]].EdgeID.push_back(EDGES[i].id);
		G[hashTableVertex[EDGES[i].node2id]].AdjList.push_back(hashTableVertex[EDGES[i].node1id]);
		G[hashTableVertex[EDGES[i].node2id]].EdgeID.push_back(EDGES[i].id);
	}

	// CREATE CYCLE LIST AND FIND ALL CYCLES
	vector<vector<int> > CYCLES;
	detectCycles(G, CYCLES);

	std::cout << "NUM TOTAL CYCLES: " << CYCLES.size() << std::endl;
	for (int j = 0; j < CYCLES.size(); ++j) {
		cout << "Cycle #" << j << " ";
		for (int i = 0; i < CYCLES[j].size(); ++i) {
			cout << NODES[CYCLES[j][i]].id << " ";
		}
		cout << endl;
	}

	int InitialSubsetSize = 1; // DEFINES THE NUMBER OF INITIAL SUBSET OF FEASIBLE POINTS PER SUBPROBLEM TO INITIALIZE THE MASTER PROBLEM
	bool randomizeInitialObjective = false;

	// CALCULATE NODE DEMANDS
	for (int j = 0; j < LOADS.size(); ++j) {
		//cout << LOADS[j].realphase[0] << " " << LOADS[j].realphase[1] << " " << LOADS[j].realphase[2] << " " << endl;
		for (int k = 0; k < NumPhases; ++k) {
			if (LOADS[j].hasphase[k]) {
				NODES[hashTableVertex[LOADS[j].node_id]].hasphase[k] = true;
				// Calculate Total Demand at Node 
				NODES[hashTableVertex[LOADS[j].node_id]].demand[k] += LOADS[j].realphase[k];
			}
		}
	}
	// CHECK IF NODE HAS GENERATOR
	for (int j = 0; j < GENERATORS_mod.size(); ++j) {
		for (int k = 0; k < NumPhases; ++k) {
			if (GENERATORS_mod[j].hasphase[k]) {
				NODES[hashTableVertex[GENERATORS_mod[j].node_id]].hasgenerator[k] = true;
			}
		}
	}

	vector<int> NumDamagedEdges(NumScenarios + 1);

	int ColSize = 3 * EDGES.size() + GENERATORS_mod.size();

	string input_file = root_str + "config-master.xml";
	vector<dataNode<double> > LINE_CONSTRUCTION_COST;
	vector<dataNode<double> > MICROGRID_COST_mod;
	vector<dataNode<double> > MICROGRID_FIXED_COST_mod;
	vector<dataNode<bool> > IS_CRITICAL_LOAD;
	vector<dataNode<double> > MAX_MICROGRID_mod;
	vector<dataNode<double> > HARDEN_COST;
	vector<dataNode<double> > LINE_SWITCH_COST;
	vector<dataNode<string> > BUILD_SET;
	vector<dataNode<string> > HARDEN_SET;
	double LoadMet;
	double CriticalLoadRead;
	double PhaseVariation;

	cout << input_file << endl;
	readMASTERData(input_file, 
		LINE_CONSTRUCTION_COST,
		MICROGRID_COST_mod,
		MICROGRID_FIXED_COST_mod,
		IS_CRITICAL_LOAD,
		MAX_MICROGRID_mod,
		HARDEN_COST,
		LINE_SWITCH_COST,
		BUILD_SET,
		HARDEN_SET,
		LoadMet,
		CriticalLoadRead,
		PhaseVariation);
	cout << "FINISHED MASTER DATA" << endl;

	// CREATE GRAPH
	vector<graph_vertex> G_s(NODES.size());
	for (int i = 0; i < EDGES.size(); ++i) {
		for (int l = 0; l < LINE_SWITCH_COST.size(); ++l) {
			if (i == hashTableEdge[LINE_SWITCH_COST[l].id]) {
				G_s[hashTableVertex[EDGES[i].node1id]].AdjList.push_back(hashTableVertex[EDGES[i].node2id]);
				G_s[hashTableVertex[EDGES[i].node1id]].EdgeID.push_back(EDGES[i].id);
				G_s[hashTableVertex[EDGES[i].node2id]].AdjList.push_back(hashTableVertex[EDGES[i].node1id]);
				G_s[hashTableVertex[EDGES[i].node2id]].EdgeID.push_back(EDGES[i].id);
				break;
			}
		}
	}

	// CREATE CYCLE LIST AND FIND ALL CYCLES
	vector<vector<int> > CYCLES_s;
	detectCycles(G_s, CYCLES_s);

	std::cout << "NUM TOTAL UNBREAKABLE CYCLES: " << CYCLES_s.size() << std::endl;
	for (int j = 0; j < CYCLES_s.size(); ++j) {
		cout << "Cycle #" << j << " ";
		for (int i = 0; i < CYCLES_s[j].size(); ++i) {
			cout << NODES[CYCLES_s[j][i]].id << " ";
		}
		cout << endl;
	}

	// -- ADD MICROGRID GENERATION TO ALL NODES ---------------------------------
	if (addMoreGenerators) {
		for (int i = 0; i < NODES.size(); ++i) {
			string generatorID = "gA" + NODES[i].id;
			MAX_MICROGRID_mod.push_back(dataNode<double>(generatorID, 5000.0));
			MICROGRID_FIXED_COST_mod.push_back(dataNode<double>(generatorID, 500.0));
			MICROGRID_COST_mod.push_back(dataNode<double>(generatorID, 2.7));
		}
	}
	// -- MICROGRID DISCRETIZATION ---------------------------------------------
	vector<int> facilityVariablePerFacility(GENERATORS_mod.size(), 0);
	for (int j = 0; j < MAX_MICROGRID_mod.size(); ++j) {
		//cout << hashTableGenerator[MAX_MICROGRID_mod[j].id] << " " << MAX_MICROGRID_mod[j].id << " " << MAX_MICROGRID_mod[j].data << endl;
		if (GENERATORS_mod[hashTableGenerator[MAX_MICROGRID_mod[j].id]].maxrealphase[0] < IloInfinity) {
			facilityVariablePerFacility[hashTableGenerator[MAX_MICROGRID_mod[j].id]] = int(MAX_MICROGRID_mod[j].data / discreteMicrogrid);
		}
	}

	//cout << maxelement<int>(facilityVariablePerFacility) << " " << minelement<int>(facilityVariablePerFacility) << endl;

	vector<generatorData> GENERATORS(0);
	vector<dataNode<double> > MAX_MICROGRID(0);
	vector<dataNode<double> > MICROGRID_FIXED_COST(0);
	vector<dataNode<double> > MICROGRID_COST(0);
	
	int numFacility = -1;
	for (int j = 0; j < GENERATORS_mod.size(); ++j) {
		//cout << GENERATORS_mod[j].id << " "  << facilityVariablePerFacility[j] << endl;
		if (facilityVariablePerFacility[j] == 0) {
			cout << GENERATORS_mod[j].id << " "  << GENERATORS_mod[j].node_id << " " << 1 << endl;
			numFacility++;
			GENERATORS.push_back(generatorData());
			GENERATORS[numFacility].id = GENERATORS_mod[j].id;
			GENERATORS[numFacility].node_id = GENERATORS_mod[j].node_id;
			for (int k = 0; k < NumPhases; ++k) {
				GENERATORS[numFacility].hasphase.push_back(GENERATORS_mod[j].hasphase[k]);
				GENERATORS[numFacility].maxrealphase.push_back(GENERATORS_mod[j].maxrealphase[k]);
				GENERATORS[numFacility].maxreactivephase.push_back(GENERATORS_mod[j].maxreactivephase[k]);
			}
		} else {
			cout << GENERATORS_mod[j].id << " "  << GENERATORS_mod[j].node_id << " " << facilityVariablePerFacility[j] << endl;
			for (int l = 0; l < facilityVariablePerFacility[j]; ++l) {
				std::stringstream ss1;
				ss1 << l;	
				string newID = GENERATORS_mod[j].id + "_" + ss1.str();
				numFacility++;
				GENERATORS.push_back(generatorData());
				GENERATORS[numFacility].id = newID;
				GENERATORS[numFacility].node_id = GENERATORS_mod[j].node_id;
				for (int k = 0; k < NumPhases; ++k) {
					GENERATORS[numFacility].hasphase.push_back(GENERATORS_mod[j].hasphase[k]);
					GENERATORS[numFacility].maxrealphase.push_back(GENERATORS_mod[j].maxrealphase[k] / (double) facilityVariablePerFacility[j]);
					GENERATORS[numFacility].maxreactivephase.push_back(GENERATORS_mod[j].maxreactivephase[k] / (double) facilityVariablePerFacility[j]);
				}
			}
		}
	}
	
	for (int j = 0; j < MAX_MICROGRID_mod.size(); ++j) {
		for (int l = 0; l < facilityVariablePerFacility[hashTableGenerator[MAX_MICROGRID_mod[j].id]]; ++l) {
			std::stringstream ss1;
			ss1 << l;	
			string newID = MAX_MICROGRID_mod[j].id + "_" + ss1.str();
			MAX_MICROGRID.push_back(dataNode<double>(newID, (l+1) * discreteMicrogrid));
		}
	}
	for (int j = 0; j < MICROGRID_COST_mod.size(); ++j) {
		for (int l = 0; l < facilityVariablePerFacility[hashTableGenerator[MICROGRID_COST_mod[j].id]]; ++l) {
			std::stringstream ss1;
			ss1 << l;	
			string newID = MICROGRID_COST_mod[j].id + "_" + ss1.str();
			MICROGRID_COST.push_back(dataNode<double>(newID, MICROGRID_COST_mod[j].data));
		}
	}
	for (int j = 0; j < MICROGRID_FIXED_COST_mod.size(); ++j) {
		for (int l = 0; l < facilityVariablePerFacility[hashTableGenerator[MICROGRID_FIXED_COST_mod[j].id]]; ++l) {
			std::stringstream ss1;
			ss1 << l;	
			string newID = MICROGRID_FIXED_COST_mod[j].id + "_" + ss1.str();
			MICROGRID_FIXED_COST.push_back(dataNode<double>(newID, MICROGRID_FIXED_COST_mod[j].data));
		}
	}

	GENERATORS_mod.clear();
	MAX_MICROGRID_mod.clear();
	MICROGRID_COST_mod.clear();
	MICROGRID_FIXED_COST_mod.clear();

	hashTableGenerator.clear();
	for (int i = 0; i < GENERATORS.size(); ++i) {
		hashTableGenerator.insert(pair<string, int>(GENERATORS[i].id, i));
	}

	ColSize = 3 * EDGES.size() + GENERATORS.size();

	cout << "NUM CRITICAL LOADS: " << IS_CRITICAL_LOAD.size() << endl;
	// SCENARIO DATA
	vector<vector<dataNode<bool> > > DISABLED(NumScenarios+1,vector<dataNode<bool> >());
	vector<vector<dataNode<bool> > > HARDENED_DISABLED(NumScenarios+1, vector<dataNode<bool> >());
	// FIX THE SEED
	cout << "DAMAGE CALCULATION STARTED..." << endl;
	srand(0);
	double damage_probability = atof(damage_str.c_str());
	// READ SCENARIO DATA
	for (int i = 0; i < NumScenarios + 1; ++i) {
		int numDamage = 0;
		for (int j = 0; j < EDGES.size(); ++j) {
			if (i != NumScenarios) {
				bool dmg_ind = false;
				for (int k = 0; k < EDGES[j].NumPoles; ++k) {
					double dice_roll = (double) rand() / RAND_MAX;
					//cout << "RANDOM ROLL " << dice_roll << " " << damage_probability << endl;
					if (dice_roll < damage_probability) {
						dmg_ind = true;
						numDamage++;
						break;
					}
				}
				DISABLED[i].push_back(dataNode<bool>(EDGES[j].id, dmg_ind));	
				HARDENED_DISABLED[i].push_back(dataNode<bool>(EDGES[j].id, false));	
			} else {
				DISABLED[i].push_back(dataNode<bool>(EDGES[j].id, false));	
				HARDENED_DISABLED[i].push_back(dataNode<bool>(EDGES[j].id, false));	
			}
		}
		cout << "Scenario " << i << " numDamage " << numDamage << endl;
	}

	cout << "DAMAGE CREATED" << endl;

	vector<double> ObjectiveVector(ColSize,0.0);

	for (int i = 0; i < NumScenarios + 1; ++i) {

		PopulateScenarioData(
			SubSolver[i],
			hashTableVertex,
			hashTableEdge,
			hashTableGenerator,
			hashTableLoad,
			NODES,
			EDGES,
			GENERATORS,
			LOADS,
			LINECODES,
			LINE_CONSTRUCTION_COST,
			MICROGRID_COST,
			MICROGRID_FIXED_COST,
			IS_CRITICAL_LOAD,
			MAX_MICROGRID,
			HARDEN_COST,
			LINE_SWITCH_COST,
			LoadMet,
			CriticalLoadMet,
			PhaseVariation,
			G,
			CYCLES,
			SubObjective[i],
			ObjectiveVector,
			lineExistsConstraint[i], // CONSTRAINTS
			lineConstraint1[i], // CONSTRAINTS
			lineConstraint2[i],
			lineReactiveConstraint1[i],
			lineReactiveConstraint2[i],
			directionConstraint[i],
			switchConstraint1[i],
			switchConstraint2[i],
			switchReactiveConstraint1[i],
			switchReactiveConstraint2[i],
			variationConstraint1[i],
			variationConstraint2[i],
			variationReactiveConstraint1[i],
			variationReactiveConstraint2[i],
			damageConstraint[i],
			loadConstraint[i],
			loadReactiveConstraint[i],
			generationConstraint[i],
			generationReactiveConstraint[i],
			balanceConstraint[i],
			balanceReactiveConstraint[i],
			microgridConstraint[i],
			cycleConstraint[i],
			cycleSubConstraint[i],
			switchableConstraint[i],
			criticalServeConstraint[i],
			criticalServeReactiveConstraint[i],
			loadServeConstraint[i],
			loadServeReactiveConstraint[i],
			switchImplicationConstraint1[i],
			//switchImplicationConstraint2[i],
			hardenConstraint[i],
			voltageConstraint[i],
			voltageOffsetConstraint[i],
			voltageVariationConstraint1[i],
			voltageVariationConstraint2[i],
			lineUseVariable[i], // VARIABLES
			lineExistsVariable[i], // VARIABLES
			lineDirectionVariable[i],
			lineCycleVariable[i],
			lineHardenVariable[i],
			switchUseVariable[i],
			switchCycleVariable[i],
			flowVariable[i],
			flowReactiveVariable[i],
			voltageVariable[i],
			voltageOffsetVariable[i],
			loadVariable[i],
			loadReactiveVariable[i],
			loadServeVariable[i],
			generatorVariable[i],
			generatorReactiveVariable[i],
			//microgridVariable[i],
			facilityVariable[i],
			NumDamagedEdges[i],
			NumPhases,
			i,
			NumScenarios,
			discreteMicrogrid,
			damage_str,
			root_str,
			HardenedDisabled,
			DISABLED[i],
			HARDENED_DISABLED[i],
			newImpedanceMultiplier,
			LDFIndicator);

		// -- ADD ALL OF THE ABOVE OBJECTS MODEL -----------------------------------------------------------------
		SubProblem[i].add(lineExistsConstraint[i]);
		SubProblem[i].add(lineConstraint1[i]);
		SubProblem[i].add(lineConstraint2[i]);
		SubProblem[i].add(lineReactiveConstraint1[i]);
		SubProblem[i].add(lineReactiveConstraint2[i]);
		SubProblem[i].add(directionConstraint[i]);
		SubProblem[i].add(switchConstraint1[i]);
		SubProblem[i].add(switchConstraint2[i]);
		SubProblem[i].add(switchReactiveConstraint1[i]);
		SubProblem[i].add(switchReactiveConstraint2[i]);
		SubProblem[i].add(variationConstraint1[i]);
		SubProblem[i].add(variationConstraint2[i]);
		SubProblem[i].add(variationReactiveConstraint1[i]);
		SubProblem[i].add(variationReactiveConstraint2[i]);
		SubProblem[i].add(damageConstraint[i]);
		SubProblem[i].add(loadConstraint[i]);
		SubProblem[i].add(loadReactiveConstraint[i]);
		SubProblem[i].add(generationConstraint[i]);
		SubProblem[i].add(generationReactiveConstraint[i]);
		SubProblem[i].add(balanceConstraint[i]);
		SubProblem[i].add(balanceReactiveConstraint[i]);
		SubProblem[i].add(microgridConstraint[i]);
		SubProblem[i].add(cycleConstraint[i]);
		SubProblem[i].add(cycleSubConstraint[i]);
		SubProblem[i].add(switchableConstraint[i]);
		SubProblem[i].add(criticalServeConstraint[i]);
		SubProblem[i].add(criticalServeReactiveConstraint[i]);
		SubProblem[i].add(loadServeConstraint[i]);
		SubProblem[i].add(loadServeReactiveConstraint[i]);
		SubProblem[i].add(switchImplicationConstraint1[i]);
		//SubProblem[i].add(switchImplicationConstraint2[i]);
		SubProblem[i].add(hardenConstraint[i]);
		SubProblem[i].add(voltageConstraint[i]);
		SubProblem[i].add(voltageOffsetConstraint[i]);
		SubProblem[i].add(voltageVariationConstraint1[i]);
		SubProblem[i].add(voltageVariationConstraint2[i]);

		SubProblem[i].add(lineUseVariable[i]);
		SubProblem[i].add(lineExistsVariable[i]);
		SubProblem[i].add(lineDirectionVariable[i]);
		SubProblem[i].add(lineCycleVariable[i]);
		SubProblem[i].add(lineHardenVariable[i]);
		SubProblem[i].add(switchUseVariable[i]);
		SubProblem[i].add(switchCycleVariable[i]);
		SubProblem[i].add(flowVariable[i]);
		SubProblem[i].add(flowReactiveVariable[i]);
		SubProblem[i].add(voltageVariable[i]);
		SubProblem[i].add(voltageOffsetVariable[i]);
		SubProblem[i].add(loadVariable[i]);
		SubProblem[i].add(loadReactiveVariable[i]);
		SubProblem[i].add(loadServeVariable[i]);
		SubProblem[i].add(generatorVariable[i]);
		SubProblem[i].add(generatorReactiveVariable[i]);
		//SubProblem[i].add(microgridVariable[i]);
		SubProblem[i].add(facilityVariable[i]);

		SubProblem[i].add(SubObjective[i]);
		/*for (int j = 0; j <= i; ++j) { 
			cout << "OBJECTIVE SIZE: " << SubObjective[j].getSize() << endl;
		}*/
		// -- EXTRACT THE PROBLEM --------------------------------------------------------------------------------------
		SubSolver[i].cplex.extract(SubProblem[i]);

		//SubSolver[i].cplex.exportModel("subproblem.lp");
		//double MAX_OBJ_COEF = 1e+3;

		if (SOLVE) {

			bool outputResults = false;

			// -- OBJECTIVE -------------------------------------------------------------------------------------------

			SubSolver[i].cplex.solve();
			cout << "SOLVED! " << SubSolver[i].cplex.getStatus() << " OBJECTIVE: " << SubSolver[i].cplex.getObjValue() << endl;

			IloEnv newEnv = SubSolver[i].cplex.getEnv();
			IloNumArray voltagevals(newEnv);
			IloNumArray lineusevals(newEnv);
			IloNumArray switchusevals(newEnv);
			IloNumArray linehardenvals(newEnv);
			IloNumArray voltageoffsetvals(newEnv);
			IloNumArray flowvals(newEnv);
			IloNumArray flowreactvals(newEnv);
			IloNumArray loadvals(newEnv);
			SubSolver[i].cplex.getValues(voltagevals, voltageVariable[i]);
			SubSolver[i].cplex.getValues(lineusevals, lineUseVariable[i]);
			SubSolver[i].cplex.getValues(switchusevals, switchUseVariable[i]);
			SubSolver[i].cplex.getValues(linehardenvals, lineHardenVariable[i]);
			SubSolver[i].cplex.getValues(voltageoffsetvals, voltageOffsetVariable[i]);
			SubSolver[i].cplex.getValues(flowvals, flowVariable[i]);
			SubSolver[i].cplex.getValues(flowreactvals, flowReactiveVariable[i]);
			SubSolver[i].cplex.getValues(loadvals, loadVariable[i]);

			for (int j = 0; j < NODES.size(); ++j) {
				cout << NODES[j].id;
				for (int k = 0; k < NumPhases; ++k) {
					cout << " " << (NODES[j].hasphase[k] ? "true" : "false");
				}
				cout << " (";
				for (int k = 0; k < NumPhases; ++k) {
					cout << " " << pow(voltagevals[j * NumPhases + k] / 1e+3, 0.5);
				}
				for (int k = 0; k < NumPhases; ++k) {
					cout << " " << loadvals[j * NumPhases + k];
				}

				cout << " ) " << endl;
			}
			for (int j = 0; j < EDGES.size(); ++j) {
				cout << EDGES[j].id << " " << EDGES[j].node1id << " "  << EDGES[j].node2id << " " << EDGES[j].length << " (" << lineusevals[j] << " " << switchusevals[j] << " " << linehardenvals[j];
				for (int k = 0; k < NumPhases; ++k) {
					cout << " " << voltageoffsetvals[j * NumPhases + k];
				}
				for (int k = 0; k < NumPhases; ++k) {
					cout << " " << round(flowvals[j * NumPhases + k]);
				}
				for (int k = 0; k < NumPhases; ++k) {
					cout << " " << round(flowreactvals[j * NumPhases + k]);
				}
				cout << " ) " << endl;
			}
			

			if (SubSolver[i].cplex.getStatus() == IloAlgorithm::Infeasible) {
			// -- CHECK CONFLICTS --------------------------------------------------------------------------------------

				IloConstraintArray check(env);
				check.add(lineExistsConstraint[i]);
				check.add(lineConstraint1[i]);
				check.add(lineConstraint2[i]);
				check.add(lineReactiveConstraint1[i]);
				check.add(lineReactiveConstraint2[i]);
				check.add(directionConstraint[i]);
				check.add(switchConstraint1[i]);
				check.add(switchConstraint2[i]);
				check.add(switchReactiveConstraint1[i]);
				check.add(switchReactiveConstraint2[i]);
				check.add(variationConstraint1[i]);
				check.add(variationConstraint2[i]);
				check.add(variationReactiveConstraint1[i]);
				check.add(variationReactiveConstraint2[i]);
				check.add(damageConstraint[i]);
				check.add(loadConstraint[i]);
				check.add(loadReactiveConstraint[i]);
				check.add(generationConstraint[i]);
				check.add(generationReactiveConstraint[i]);
				check.add(balanceConstraint[i]);
				check.add(balanceReactiveConstraint[i]);
				check.add(microgridConstraint[i]);
				check.add(cycleConstraint[i]);
				check.add(cycleSubConstraint[i]);
				check.add(switchableConstraint[i]);
				check.add(criticalServeConstraint[i]);
				check.add(criticalServeReactiveConstraint[i]);
				check.add(loadServeConstraint[i]);
				check.add(loadServeReactiveConstraint[i]);
				check.add(switchImplicationConstraint1[i]);
				//check.add(switchImplicationConstraint2[i]);
				check.add(hardenConstraint[i]);
				check.add(voltageConstraint[i]);
				check.add(voltageOffsetConstraint[i]);
				check.add(voltageVariationConstraint1[i]);
				check.add(voltageVariationConstraint2[i]);

				IloNumArray prefs(env);
				for (int j = 0; j < check.getSize(); ++j) {
					prefs.add(1.0);
				}
				if (SubSolver[i].cplex.refineConflict(check, prefs)) {
					std::ofstream conflictOutput("conflicts.txt");
					IloCplex::ConflictStatusArray conflict = SubSolver[i].cplex.getConflict(check);
					env.getImpl()->useDetailedDisplay(IloTrue);
					conflictOutput << "Conflict: " << endl;
					for (int j = 0; j < check.getSize(); ++j) {
						if (conflict[j] == IloCplex::ConflictMember) {
							conflictOutput << "Proved: " << check[j] << endl;
						} 
						if (conflict[j] == IloCplex::ConflictPossibleMember) {
							conflictOutput << "Possible: " << check[j] << endl;
						}
					}
				}
			} else {
				env = SubSolver[i].cplex.getEnv();
				IloNumArray lineusevals(env);
				IloNumArray lineswitchvals(env);
				IloNumArray linehardenvals(env);
				IloNumArray facilityvals(env);

				SubSolver[i].cplex.getValues(lineusevals, lineUseVariable[i]);
				SubSolver[i].cplex.getValues(lineswitchvals, switchUseVariable[i]);
				SubSolver[i].cplex.getValues(linehardenvals, lineHardenVariable[i]);
				//SubSolver[i].cplex.getValues(microgridvals, microgridVariable[i]);
				SubSolver[i].cplex.getValues(facilityvals, facilityVariable[i]);
				// cout << "DID IT!!" << endl;
				// int ColSize = 3 * EDGES.size() + (NumPhases + 1) * GENERATORS.size();
				// vector<vector<vector<double> > > SubColumn (NumScenarios, vector<vector<double> >(0));
				SubColumn[i].push_back(vector<double>(ColSize));
				for (int j = 0; j < EDGES.size(); ++j) {
					SubColumn[i][SubColumn[i].size()-1][j] = lineusevals[j];
					SubColumn[i][SubColumn[i].size()-1][EDGES.size() + j] = lineswitchvals[j];
					SubColumn[i][SubColumn[i].size()-1][2*EDGES.size() + j] = linehardenvals[j];
				}
				for (int j = 0; j < GENERATORS.size(); ++j) {
					SubColumn[i][SubColumn[i].size()-1][3 * EDGES.size() + j] = facilityvals[j];
				}

				if (outputResults) { // Output solution
					//cout << "SOLUTION STATUS: " << cplex.getCplexStatus() << endl;
					//cout << "OBJECTIVE VALUE: " << cplex.getObjValue() << endl;

					IloNumArray flowvals(env);
					IloNumArray lineusevals(env);
					IloNumArray lineswitchvals(env);
					SubSolver[i].cplex.getValues(lineusevals, lineUseVariable[i]);
					SubSolver[i].cplex.getValues(flowvals, flowVariable[i]);
					SubSolver[i].cplex.getValues(lineswitchvals, switchUseVariable[i]);
					cout << "FLOW VALUES: " << endl;
					for (int j = 0; j < EDGES.size(); ++j) {
						//cout << "EDGE " << EDGES[j].id << ": " << lineusevals[j] << endl; 
						cout << "EDGE " << EDGES[j].id << " (" << lineusevals[j] << "," << lineswitchvals[j] << "," << linehardenvals[j] << ") VALUES: (" << flowvals[j*NumPhases]  << "," << flowvals[j*NumPhases+1]  << "," << flowvals[j*NumPhases+2]  << ")" << endl;
					}
					IloNumArray genvals(env);
					IloNumArray facvals(env);
					SubSolver[i].cplex.getValues(genvals, generatorVariable[i]);
					SubSolver[i].cplex.getValues(facvals, facilityVariable[i]);
					cout << "GENERATION VALUES: " << endl;
					for (int j = 0; j < GENERATORS.size(); ++j) {
						double sum = 0.0;
						for (int k = 0; k < NumPhases; ++k) {
							sum += genvals[j*NumPhases + k];
						}
						cout << "GENERATOR #" << GENERATORS[j].id << " (" << facvals[j] << "," << sum << ")" << endl;  
					}

					/*for (int j = 0; j < EDGES.size(); ++j) {
						if (lineusevals[j] == 1) {
							cout << EDGES[j].node1id << " " << EDGES[j].node2id << endl;
						}
					}*/
	
					cout << "SOLUTION STATUS: " << SubSolver[i].cplex.getCplexStatus() << endl;
					cout << "OBJECTIVE VALUE: " << SubSolver[i].cplex.getObjValue() << endl;
					/*try {
						SubSolver[i].cplex.exportModel("lp1.lp");
					}
					catch (IloException &e) {
						cerr << "Concert exception caught: " << e << endl;
					} catch (...) {
						cerr << "Unknown exception caught: " << endl;
					}*/

					// WRITE GRAPH FILES
					if (i == NumScenarios) {
					string graph_node_output = "graph/nodes_" + extent_str + ".txt"; 
					std::ofstream outputNodes(graph_node_output.c_str());
					//outputGraph << NODES.size() << " " << EDGES.size() << std::endl;
					for (int j = 0; j < NODES.size(); ++j) {
						outputNodes << NODES[j].x << " " << NODES[j].y << " " << NODES[j].hasgenerator[0] << " " << NODES[j].hasgenerator[1]<< " " << NODES[j].hasgenerator[2] << std::endl;
					}
					outputNodes.close();
					string graph_edge_output = "graph/edges_" + extent_str + ".txt";
					std::ofstream outputEdges(graph_edge_output.c_str());
					for (int j = 0; j < EDGES.size(); ++j) {
						int existing = (EDGES[j].id[0] == 'n' || EDGES[j].id[0] == 'o') ? 0 : 1;
						outputEdges << hashTableVertex[EDGES[j].node1id] << " " << hashTableVertex[EDGES[j].node2id] << " " << existing << std::endl;
					}
					outputEdges.close();
					}
		
					// CHECK IF THERE ARE CYCLES IN THE SOLUTION:
					bool CHECK_CYCLES = true;
					if (CHECK_CYCLES) {
						vector<bool> edgeIncluded (EDGES.size(), false);
						for (int j = 0; j < EDGES.size(); ++j) {
							if (lineusevals[j] == 1 && lineswitchvals[j] != 1) {
								edgeIncluded[j] = true;
							}
						}
						// CREATE GRAPH
						vector<graph_vertex> GSub(NODES.size());
						for (int j = 0; j < EDGES.size(); ++j) {
							if (edgeIncluded[j]) {
								GSub[hashTableVertex[EDGES[j].node1id]].AdjList.push_back(hashTableVertex[EDGES[j].node2id]);
								GSub[hashTableVertex[EDGES[j].node1id]].EdgeID.push_back(EDGES[j].id);
								GSub[hashTableVertex[EDGES[j].node2id]].AdjList.push_back(hashTableVertex[EDGES[j].node1id]);
								GSub[hashTableVertex[EDGES[j].node2id]].EdgeID.push_back(EDGES[j].id);
							}
						}

						// CREATE CYCLE LIST AND FIND ALL CYCLES
						vector<vector<int> > CYCLESSub;
						detectCycles(GSub, CYCLESSub);
						cout << "SUBGRAPH HAS " << CYCLESSub.size() << " CYCLES" << endl;
					}
				}
			}
		}

		cout << "Scenario optimization ended." << endl;

		//cplex.clearModel();
		//cplex.clearCuts();
	}

	/*for (int j = 0; j < GENERATORS.size(); ++j) {
		cout << MAX_ADDED[j] << " ";
	}
	cout << endl;*/

	// COMPUTE MINIMAL FEASIBLE VECTOR:
	vector<double> MIN_FEASIBLE_VECTOR (ColSize, 0.0);
	if (SOLVE) {	
		for (int j = 0; j < ColSize; ++j) {
			double maxValue = 0.0;
			for (int i = 0; i < NumScenarios; ++i) {
				for (int k = 0; k < SubColumn[i].size(); ++k) {
					if (SubColumn[i][k][j] > maxValue) {
						maxValue = SubColumn[i][k][j];
					}
				}
			}
			MIN_FEASIBLE_VECTOR[j] = maxValue;
		}
		for (int j = 0; j < EDGES.size(); ++j) {
			double switchInUse = 0.0;
			for (int i = 0; i < NumScenarios; ++i) {
				if (MIN_FEASIBLE_VECTOR[j] - SubColumn[i][0][j] > 0.5) {
					switchInUse = 1.0;
					break;
				}
			}
			MIN_FEASIBLE_VECTOR[j+EDGES.size()] = switchInUse;
		}
		cout << "MINIMAL FEASIBLE SOLUTION OBJECTIVE: " << std::inner_product(MIN_FEASIBLE_VECTOR.begin(), MIN_FEASIBLE_VECTOR.end(), ObjectiveVector.begin(), 0.0) << endl;
	}
	double MIN_FEAS_OBJ = std::inner_product(MIN_FEASIBLE_VECTOR.begin(), MIN_FEASIBLE_VECTOR.end(), ObjectiveVector.begin(), 0.0);
	
	cout << "BENDERS STARTED" << endl;

	BendersFullScenarioGeneration(
		ObjectiveVector,
		MIN_FEAS_OBJ,
		MIN_FEASIBLE_VECTOR,
		SubProblem,
		SubSolver,
		hashTableVertex,
		hashTableEdge,
		hashTableGenerator,
		hashTableLoad,
		NODES,
		EDGES,
		GENERATORS,
		LOADS,
		LINECODES,
		LINE_CONSTRUCTION_COST,
		MICROGRID_COST,
		MICROGRID_FIXED_COST,
		IS_CRITICAL_LOAD,
		MAX_MICROGRID,
		HARDEN_COST,
		LINE_SWITCH_COST,
		LoadMet,
		CriticalLoadMet,
		PhaseVariation,
		G,
		CYCLES,
		SubObjective, // OBJECTIVE 
		lineExistsConstraint, // CONSTRAINTS
		lineConstraint1, // CONSTRAINTS
		lineConstraint2,
		lineReactiveConstraint1,
		lineReactiveConstraint2,
		directionConstraint,
		switchConstraint1,
		switchConstraint2,
		switchReactiveConstraint1,
		switchReactiveConstraint2,
		variationConstraint1,
		variationConstraint2,
		variationReactiveConstraint1,
		variationReactiveConstraint2,
		damageConstraint,
		loadConstraint,
		loadReactiveConstraint,
		generationConstraint,
		generationReactiveConstraint,
		balanceConstraint,
		balanceReactiveConstraint,
		microgridConstraint,
		cycleConstraint,
		cycleSubConstraint,
		switchableConstraint,
		criticalServeConstraint,
		criticalServeReactiveConstraint,
		loadServeConstraint,
		loadServeReactiveConstraint,
		switchImplicationConstraint1,
		//switchImplicationConstraint2,
		hardenConstraint,
		voltageConstraint,
		voltageOffsetConstraint,
		voltageVariationConstraint1,
		voltageVariationConstraint2,
		lineUseVariable, // VARIABLES
		lineExistsVariable, // VARIABLES
		lineDirectionVariable,
		lineCycleVariable,
		lineHardenVariable,
		switchUseVariable,
		switchCycleVariable,
		flowVariable,
		flowReactiveVariable,
		voltageVariable,
		voltageOffsetVariable,
		loadVariable,
		loadReactiveVariable,
		loadServeVariable,
		generatorVariable,
		generatorReactiveVariable,
		//microgridVariable,
		facilityVariable,
		NumDamagedEdges,
		NumPhases,
		EDGES.size(),
		LOADS.size(),
		GENERATORS.size(),
		NumScenarios,
		discreteMicrogrid,
		damage_str,
		extent_str,
		root_str,
		HardenedDisabled,
		ExactSolver,
		LDFIndicator,
		MaxVoltageRounds,
		CHECK_RESULTS,
		DISABLED,
		HARDENED_DISABLED,
		newImpedanceMultiplier);
	return;
	
}

// MAIN PROGRAM
int main(int argc, char **argv) {

	string damage_str;
	string extent_str;
	string root_str;
	double CriticalLoadMet;
	double discreteMicrogrid;
	bool HardenedDisabled;
	bool ExactSolver;
	bool LDFIndicator;
	int MaxVoltageRounds;
	bool CHECK_RESULTS;
	double newImpedanceMultiplier;
	int NumScenarios; // Number of Disasters

	if (argc > 1) {
		damage_str = std::string(argv[1]);
		NumScenarios = atoi(argv[2]);
		newImpedanceMultiplier = atof(argv[3]);
		CriticalLoadMet = atof(argv[4]);
		discreteMicrogrid = atof(argv[5]);
		extent_str = std::string(argv[6]);
		root_str = std::string(argv[7]) + extent_str + "/";
		HardenedDisabled = (std::string(argv[8]).compare("HardenedDisabled") == 0);
		ExactSolver = (std::string(argv[9]).compare("Exact") == 0);
		LDFIndicator = (std::string(argv[10]).compare("LDF") == 0);
		MaxVoltageRounds = atoi(argv[11]);
		CHECK_RESULTS = (std::string(argv[12]).compare("Output") == 0);
	} else {
		damage_str = "1.0% Per Mile Damage";
		NumScenarios = 100;
		newImpedanceMultiplier = 1.0;
		CriticalLoadMet = 0.98;
		discreteMicrogrid = 500.0;
		extent_str = "Rural";
		root_str = "LargeExpansion/Rural/";
		HardenedDisabled = false;
		ExactSolver = false;
		LDFIndicator = false;
		MaxVoltageRounds = 1;
		CHECK_RESULTS = false;
	}

	cout << damage_str << endl;
	cout << "IMPEDANCE MULTIPLIER: " << newImpedanceMultiplier << endl;
	cout << "CRITICAL LOAD MET: " << CriticalLoadMet << endl;
	cout << "MICROGRID DISCRETIZATION: " << discreteMicrogrid << endl;
	cout << "HARDENED DISABLED: " << (HardenedDisabled ? "true" : "false") << endl;
	cout << "EXACT SOLVER: " << (ExactSolver ? "true" : "false") << endl;
	cout << "FLOW TYPE: " << (LDFIndicator ? "LDF" : "MCF") << endl;
	cout << "OUTPUT RESULTS: " << (CHECK_RESULTS ? "true" : "false") << endl;
	cout << "MAX VOLTAGE CUTS: " << MaxVoltageRounds << endl;
	//IloEnv env;
	try {
		int SubNumCon = 0;
		
		// SUB PROBLEMS
		vector<IloModel> SubProblem (NumScenarios + 1);
		
		// Defined one by one to enhance visibility of model, we can reduce below into a single vector.
		// The names are taken from the paper for ease of use.
		vector<IloObjective> SubObjective (NumScenarios);
		
		vector<IloNumVarArray> lineUseVariable (NumScenarios + 1);
		vector<IloNumVarArray> lineExistsVariable (NumScenarios + 1);
		vector<IloNumVarArray> lineDirectionVariable (NumScenarios + 1);
		vector<IloNumVarArray> lineCycleVariable (NumScenarios + 1);
		vector<IloNumVarArray> switchCycleVariable (NumScenarios + 1);
		vector<IloNumVarArray> lineHardenVariable (NumScenarios + 1);		
		vector<IloNumVarArray> switchUseVariable (NumScenarios + 1);
		vector<IloNumVarArray> flowVariable (NumScenarios + 1);
		vector<IloNumVarArray> flowReactiveVariable (NumScenarios + 1);
		vector<IloNumVarArray> voltageVariable (NumScenarios + 1);
		vector<IloNumVarArray> voltageOffsetVariable (NumScenarios + 1);
		vector<IloNumVarArray> loadVariable (NumScenarios + 1);
		vector<IloNumVarArray> loadReactiveVariable (NumScenarios + 1);
		vector<IloNumVarArray> loadServeVariable (NumScenarios + 1);
		vector<IloNumVarArray> generatorVariable (NumScenarios + 1);		
		vector<IloNumVarArray> generatorReactiveVariable (NumScenarios + 1);		
		//vector<IloNumVarArray> microgridVariable (NumScenarios + 1);
		vector<IloNumVarArray> facilityVariable (NumScenarios + 1);


		vector<IloRangeArray> lineExistsConstraint (NumScenarios + 1);
		vector<IloRangeArray> lineConstraint1 (NumScenarios + 1);
		vector<IloRangeArray> lineConstraint2 (NumScenarios + 1);
		vector<IloRangeArray> lineReactiveConstraint1 (NumScenarios + 1);
		vector<IloRangeArray> lineReactiveConstraint2 (NumScenarios + 1);
		vector<IloRangeArray> directionConstraint (NumScenarios + 1);
		vector<IloRangeArray> switchConstraint1 (NumScenarios + 1);
		vector<IloRangeArray> switchConstraint2 (NumScenarios + 1);
		vector<IloRangeArray> switchReactiveConstraint1 (NumScenarios + 1);
		vector<IloRangeArray> switchReactiveConstraint2 (NumScenarios + 1);
		vector<IloRangeArray> variationConstraint1 (NumScenarios + 1);
		vector<IloRangeArray> variationConstraint2 (NumScenarios + 1);
		vector<IloRangeArray> variationReactiveConstraint1 (NumScenarios + 1);
		vector<IloRangeArray> variationReactiveConstraint2 (NumScenarios + 1);
		vector<IloRangeArray> damageConstraint (NumScenarios + 1);		
		vector<IloRangeArray> loadConstraint (NumScenarios + 1);
		vector<IloRangeArray> loadReactiveConstraint (NumScenarios + 1);
		vector<IloRangeArray> generationConstraint (NumScenarios + 1);
		vector<IloRangeArray> generationReactiveConstraint (NumScenarios + 1);
		vector<IloRangeArray> balanceConstraint (NumScenarios + 1);		
		vector<IloRangeArray> balanceReactiveConstraint (NumScenarios + 1);		
		vector<IloRangeArray> microgridConstraint (NumScenarios + 1);
		vector<IloRangeArray> cycleConstraint (NumScenarios + 1);
		vector<IloRangeArray> cycleSubConstraint (NumScenarios + 1);
		vector<IloRangeArray> switchableConstraint (NumScenarios + 1);		
		vector<IloRangeArray> criticalServeConstraint (NumScenarios + 1);
		vector<IloRangeArray> criticalServeReactiveConstraint (NumScenarios + 1);
		vector<IloRangeArray> loadServeConstraint (NumScenarios + 1);
		vector<IloRangeArray> loadServeReactiveConstraint (NumScenarios + 1);
		vector<IloRangeArray> switchImplicationConstraint1 (NumScenarios + 1);
		//vector<IloRangeArray> switchImplicationConstraint2 (NumScenarios + 1);
		vector<IloRangeArray> hardenConstraint (NumScenarios + 1);
		vector<IloRangeArray> voltageConstraint (NumScenarios + 1);
		vector<IloRangeArray> voltageOffsetConstraint (NumScenarios + 1);
		vector<IloRangeArray> voltageVariationConstraint1 (NumScenarios + 1);
		vector<IloRangeArray> voltageVariationConstraint2 (NumScenarios + 1);


		// GET SOLVERS
      		Solver MasterSolver = Solver();
		vector<Solver> SubSolver(0);
		for (int i = 0; i < NumScenarios + 1; ++i) {
			SubSolver.push_back(Solver());
		}

		// MASTER PROBLEM // Constraints will complain!
		IloEnv env = MasterSolver.cplex.getEnv();
		IloModel MasterProblem (env);
		IloNumArray mMaster(env);
		//IloNumArray mSub(env);

		IloObjective   MasterObjective = IloAdd(MasterProblem, IloMinimize(env));
		IloRangeArray  MasterConstraints = IloAdd(MasterProblem, IloRangeArray(env, mMaster, IloInfinity));
		IloNumVarArray MasterVariables(env);

		for (int i = 0; i < NumScenarios + 1; ++i) {
			env = SubSolver[i].cplex.getEnv();
			SubProblem[i] = IloModel(env);

			SubObjective[i] = IloObjective(env);
			SubObjective[i].setSense(IloObjective::Minimize);

			lineUseVariable[i] = IloNumVarArray(env);
			lineExistsVariable[i] = IloNumVarArray(env);
			lineDirectionVariable[i] = IloNumVarArray(env);
			lineCycleVariable[i] = IloNumVarArray(env);
			switchCycleVariable[i] = IloNumVarArray(env);
			lineHardenVariable[i] = IloNumVarArray(env);
			switchUseVariable[i] = IloNumVarArray(env);
			flowVariable[i] = IloNumVarArray(env);
			flowReactiveVariable[i] = IloNumVarArray(env);
			voltageVariable[i] = IloNumVarArray(env);
			voltageOffsetVariable[i] = IloNumVarArray(env);
			loadVariable[i] = IloNumVarArray(env);
			loadReactiveVariable[i] = IloNumVarArray(env);
			loadServeVariable[i] = IloNumVarArray(env);
			generatorVariable[i] = IloNumVarArray(env);
			generatorReactiveVariable[i] = IloNumVarArray(env);
			//microgridVariable[i] = IloNumVarArray(env);
			facilityVariable[i] = IloNumVarArray(env);


			lineExistsConstraint[i] = IloRangeArray(env);
			lineConstraint1[i] = IloRangeArray(env);
			lineConstraint2[i] = IloRangeArray(env);
			lineReactiveConstraint1[i] = IloRangeArray(env);
			lineReactiveConstraint2[i] = IloRangeArray(env);
			directionConstraint[i] = IloRangeArray(env);
			switchConstraint1[i] = IloRangeArray(env);
			switchConstraint2[i] = IloRangeArray(env);
			switchReactiveConstraint1[i] = IloRangeArray(env);
			switchReactiveConstraint2[i] = IloRangeArray(env);
			variationConstraint1[i] = IloRangeArray(env);
			variationConstraint2[i] = IloRangeArray(env);
			variationReactiveConstraint1[i] = IloRangeArray(env);
			variationReactiveConstraint2[i] = IloRangeArray(env);
			damageConstraint[i] = IloRangeArray(env);		
			loadConstraint[i] = IloRangeArray(env);
			loadReactiveConstraint[i] = IloRangeArray(env);
			generationConstraint[i] = IloRangeArray(env);
			generationReactiveConstraint[i] = IloRangeArray(env);
			balanceConstraint[i] = IloRangeArray(env);		
			balanceReactiveConstraint[i] = IloRangeArray(env);		
			microgridConstraint[i] = IloRangeArray(env);
			cycleConstraint[i] = IloRangeArray(env);
			cycleSubConstraint[i] = IloRangeArray(env);
			switchableConstraint[i] = IloRangeArray(env);		
			criticalServeConstraint[i] = IloRangeArray(env);
			criticalServeReactiveConstraint[i] = IloRangeArray(env);
			loadServeConstraint[i] = IloRangeArray(env);
			loadServeReactiveConstraint[i] = IloRangeArray(env);
			switchImplicationConstraint1[i] = IloRangeArray(env);
			//switchImplicationConstraint2[i] = IloRangeArray(env);
			hardenConstraint[i] = IloRangeArray(env);
			voltageConstraint[i] = IloRangeArray(env);
			voltageOffsetConstraint[i] = IloRangeArray(env);
			voltageVariationConstraint1[i] = IloRangeArray(env);
			voltageVariationConstraint2[i] = IloRangeArray(env);

		}

		// READ DATA INTO ABOVE PROBLEMS
		generateInitialProblem (MasterProblem, SubProblem,
					MasterSolver, SubSolver,
					MasterObjective, MasterConstraints,
					SubObjective, // OBJECTIVE 
					lineExistsConstraint, // CONSTRAINTS
					lineConstraint1, // CONSTRAINTS
					lineConstraint2,
					lineReactiveConstraint1,
					lineReactiveConstraint2,
					directionConstraint,
					switchConstraint1,
					switchConstraint2,
					switchReactiveConstraint1,
					switchReactiveConstraint2,
					variationConstraint1,
					variationConstraint2,
					variationReactiveConstraint1,
					variationReactiveConstraint2,
					damageConstraint,
					loadConstraint,
					loadReactiveConstraint,
					generationConstraint,
					generationReactiveConstraint,
					balanceConstraint,
					balanceReactiveConstraint,
					microgridConstraint,
					cycleConstraint,
					cycleSubConstraint,
					switchableConstraint,
					criticalServeConstraint,
					criticalServeReactiveConstraint,
					loadServeConstraint,
					loadServeReactiveConstraint,
					switchImplicationConstraint1,
					//switchImplicationConstraint2,
					hardenConstraint,
					voltageConstraint,
					voltageOffsetConstraint,
					voltageVariationConstraint1,
					voltageVariationConstraint2,
					lineUseVariable, // VARIABLES
					lineExistsVariable, // VARIABLES
					lineDirectionVariable,
					lineCycleVariable,
					lineHardenVariable,
					switchUseVariable,
					switchCycleVariable,
					flowVariable,
					flowReactiveVariable,
					voltageVariable,
					voltageOffsetVariable,
					loadVariable,
					loadReactiveVariable,
					loadServeVariable,
					generatorVariable,
					generatorReactiveVariable,
					//microgridVariable,
					facilityVariable,					
					NumScenarios,
					CriticalLoadMet,
					discreteMicrogrid,
					damage_str,
					extent_str,
					root_str,
					HardenedDisabled,
					ExactSolver,
					LDFIndicator,
					MaxVoltageRounds,
					CHECK_RESULTS,
					newImpedanceMultiplier);

		
      		
	}
   	catch (IloException& ex) {
      		cerr << "Error: " << ex << endl;
   	}
   	catch (...) {
      		cerr << "Error" << endl;
   	}

   	//env.end();

   	return 0;
}
