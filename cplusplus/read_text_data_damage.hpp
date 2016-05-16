#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <map>
#include <numeric>
#include <set>

using namespace std;

struct nodeData {
public:
	nodeData() {
		hasphase.resize(3); 
		demand.resize(3);
		hasgenerator.resize(3);
		for(int i = 0; i < 3; ++i) {
			hasphase[i] = false;
			demand[i] = 0.0;
			hasgenerator[i] = false;
		}
	}
	string id;
	double x;
	double y;
	double minVoltage;
	double maxVoltage;
	double refVoltage;
	vector<bool> hasphase;
	vector<double> demand;
	vector<bool> hasgenerator;
};

struct generatorData {
public:
	generatorData() {}
	string id;
	string node_id;
	vector<bool> hasphase;
	vector<double> maxrealphase;
	vector<double> maxreactivephase;
	/*bool hasphaseA;
	bool hasphaseB;
	bool hasphaseC;*/
};

struct loadData {
public:
	loadData() {}
	string id;
	string node_id;
	vector<bool> hasphase;
	vector<double> realphase;
	vector<double> reactivephase;
	/*bool hasphaseA;
	bool hasphaseB;
	bool hasphaseC;*/
};

struct edgeData {
public:
	edgeData() {NumPhases = 0;istransformer = false;}
	string id;
	string node1id;
	string node2id;
	vector<bool> hasphase;
	/*bool hasphaseA;
	bool hasphaseB;
	bool hasphaseC;*/
	double capacity;
	double length;
	int NumPhases;
	bool istransformer;
	int linecode;
	int NumPoles;
};

struct lineCodeData {
public:
	lineCodeData() {}
	lineCodeData(int l, int n, vector<vector<double> > r, vector<vector<double> >x) {
		linecode = l;
		numphases = n;
		rmatrix.push_back(vector<vector<double> >(n, vector<double>(n)));
		xmatrix.push_back(vector<vector<double> >(n, vector<double>(n)));
		for (int i = 0; i < n; ++i) {
			//rmatrix[0].push_back(vector<double>(n));
			//xmatrix[0].push_back(vector<double>(n));
			for (int j = 0; j < n; ++j) {
				rmatrix[0][i][j] = r[i][j];
				xmatrix[0][i][j] = x[i][j];
			}
		}

		// 0 is 120 degrees rotation
		// 1 is 240 degrees rotation
		vector<vector<vector<double> > > rotation(2, 
			vector<vector<double> >(2, vector<double>(2)));
		rotation[0][0][0] = -0.5;
		rotation[0][0][1] = -0.866;
		rotation[0][1][0] = 0.866;
		rotation[0][1][1] = -0.5;
		rotation[1][0][0] = -0.5;
		rotation[1][0][1] = 0.866;
		rotation[1][1][0] = -0.866;
		rotation[1][1][1] = -0.5;		



		for (int k = 0; k < 2; ++k) {
			rmatrix.push_back(vector<vector<double> >(n, vector<double>(n)));
			xmatrix.push_back(vector<vector<double> >(n, vector<double>(n)));
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					vector<double> impedance(2);
					impedance[0] = r[i][j];
					impedance[1] = x[i][j];

					double newr = std::inner_product(rotation[k][0].begin(), rotation[k][0].end(), impedance.begin(), 0.0); 
					double newx = std::inner_product(rotation[k][1].begin(), rotation[k][1].end(), impedance.begin(), 0.0); 
					rmatrix[1+k][i][j] = newr;
					xmatrix[1+k][i][j] = newx;

				}
			}
		}


	}
	int linecode;
	int numphases;
	vector<vector<vector<double> > > rmatrix;
	vector<vector<vector<double> > > xmatrix;
};



void readTextData (string input_file,
	string transformer_file,
	string line_file,
	string linecode_file,
	vector<nodeData> &NODES,
	vector<generatorData> &GENERATORS,
	vector<loadData> &LOADS,
	vector<edgeData> &EDGES,
	map<int, lineCodeData> &LINECODES) {

	ifstream file(input_file.c_str());
	//file >> ws;
	string line;


	set<string> nodeIDs;
	while (!file.eof()) {
		getline(file, line);
		//cout << line << " " << line.length() << endl;
		if (line.compare("Nodes") == 0) {
			getline(file, line);
			while (true) {
				getline(file, line);
				//cout << line.length() << endl;
				if (!line.empty()) {
					//cout << line << endl;
					//id,x,y,minvoltage,maxvoltage
					nodeData node;
					
					istringstream ss(line);
					string s;
					getline( ss, node.id, ',' );
					if (nodeIDs.find(node.id) == nodeIDs.end()) {
						nodeIDs.insert(node.id);
					} else {
						cout << "BUGBUGBUGBUGBUBGUBGUBGUBGBUBG" << endl;
						cout << node.id << " is happening more than once!!!" << endl;
					}
					getline( ss, s, ',' );
					node.x = atof(s.c_str());
				    	getline( ss, s, ',' );
					node.y = atof(s.c_str());
					getline( ss, s, ',' );
					node.minVoltage = atof(s.c_str());
					getline( ss, s , ',');
					node.maxVoltage = atof(s.c_str());
					getline( ss, s , '\n');
					node.refVoltage = atof(s.c_str());
					NODES.push_back(node);	
				} else {
					cout << "NODES DONE " << NODES.size() << endl;
					break;
				}
			}			
		} else if (line.compare("Generators") == 0) {
			getline(file, line);
			while (true) {
				getline(file, line);
				//cout << line.length() << endl;
				if (!line.empty()) {
					//cout << line << endl;
					//id,nodeid,hasphaseA,hasphaseB,hasphaseC
					generatorData generator;
					
					istringstream ss(line);
					string s;
					getline( ss, generator.id, ',' );
					getline( ss, generator.node_id, ',' );
				    	getline( ss, s, ',' );
					generator.hasphase.push_back((s.compare("true") == 0));
					getline( ss, s, ',' );
					generator.hasphase.push_back((s.compare("true") == 0));
					getline( ss, s, ',');
					generator.hasphase.push_back((s.compare("true") == 0));
				    	getline( ss, s, ',' );
					generator.maxrealphase.push_back(atof(s.c_str()));
					getline( ss, s, ',' );
					generator.maxrealphase.push_back(atof(s.c_str()));
					getline( ss, s, ',');
					generator.maxrealphase.push_back(atof(s.c_str()));
				    	getline( ss, s, ',' );
					generator.maxreactivephase.push_back(atof(s.c_str()));
					getline( ss, s, ',' );
					generator.maxreactivephase.push_back(atof(s.c_str()));
					getline( ss, s, '\n');
					generator.maxreactivephase.push_back(atof(s.c_str()));
					GENERATORS.push_back(generator);	
				} else {
					cout << "GENERATORS DONE " << GENERATORS.size() << endl;
					break;
				}
			}			
		} else if (line.compare("Loads") == 0) {
			getline(file, line);
			while (true) {
				getline(file, line);
				//cout << line.length() << endl;
				if (!line.empty()) {
					//cout << line << endl;
					//id,nodeid,hasphaseA,hasphaseB,hasphaseC
					loadData load;
					
					istringstream ss(line);
					string s;
					getline( ss, load.id, ',' );
					getline( ss, load.node_id, ',' );
				    	getline( ss, s, ',' );
					load.hasphase.push_back((s.compare("true") == 0));
					getline( ss, s, ',' );
					load.hasphase.push_back((s.compare("true") == 0));
					getline( ss, s, ',');
					load.hasphase.push_back((s.compare("true") == 0));
				    	getline( ss, s, ',' );
					load.realphase.push_back(atof(s.c_str()));
					getline( ss, s, ',' );
					load.realphase.push_back(atof(s.c_str()));
					getline( ss, s, ',');
					load.realphase.push_back(atof(s.c_str()));
				    	getline( ss, s, ',' );
					load.reactivephase.push_back(atof(s.c_str()));
					getline( ss, s, ',' );
					load.reactivephase.push_back(atof(s.c_str()));
					getline( ss, s, '\n');
					load.reactivephase.push_back(atof(s.c_str()));
					LOADS.push_back(load);	
				} else {
					cout << "LOADS DONE " << LOADS.size() << endl;
					break;
				}
			}				
		} else if (line.compare("Edges") == 0) {
			getline(file, line);
			set<string> edgeIDs;
			while (true) {
				getline(file, line);
				//cout << line.length() << endl;
				if (!line.empty()) {
					//cout << line << endl;
					//id,nodeid,hasphaseA,hasphaseB,hasphaseC
					edgeData edge;
					
					istringstream ss(line);
					string s;
					getline( ss, edge.id, ',' );
					if (edgeIDs.find(edge.id) == edgeIDs.end()) {
						edgeIDs.insert(edge.id);
					} else {
						cout << "BUGBUGBUGBUGBUBGUBGUBGUBGBUBG" << endl;
						cout << edge.id << " is happening more than once!!!" << endl;
					}
					getline( ss, edge.node1id, ',' );
					getline( ss, edge.node2id, ',' );
					if (nodeIDs.find(edge.node1id) == nodeIDs.end() || nodeIDs.find(edge.node2id) == nodeIDs.end()) {	
						cout << "BUGBUGBUGBUGBUBGUBGUBGUBGBUBG" << endl;
						cout << edge.id << " is connected to nonexistent node!!!" << endl;
					}
				    	getline( ss, s, ',' );
					edge.hasphase.push_back((s.compare("true") == 0));
					if(edge.hasphase[0]) edge.NumPhases++;
					getline( ss, s, ',' );
					edge.hasphase.push_back((s.compare("true") == 0));
					if(edge.hasphase[1]) edge.NumPhases++;
					getline( ss, s, ',');
					edge.hasphase.push_back((s.compare("true") == 0));
					if(edge.hasphase[2]) edge.NumPhases++;
					for (int j = 0; j < 9; ++j) {
						// iterate until capacity field
						getline(ss, s, ',');
					}
					edge.capacity = atof(s.c_str());
					//edge.capacity = 15000.0;
					// get number of poles
					for (int j = 0; j < 3; ++j) {
						// iterate until numpoles field
						getline(ss, s, ',');
					}
					getline(ss, s, ',');
					edge.NumPoles = atoi(s.c_str());
					EDGES.push_back(edge);	
				} else {
					cout << "EDGES DONE " << EDGES.size() << endl;
					break;
				}
			}			
		}
	}
	file.close();

	map<string, int> hashTableEdge;
	for (int i = 0; i < EDGES.size(); ++i) {
		hashTableEdge.insert(pair<string, int>(EDGES[i].id, i));
	}

	ifstream file3(line_file.c_str());
	while (!file3.eof()) {
		getline(file3, line);
		if (line.size() > 0 && line[0] == 'N') {
			istringstream ss(line);
			string lineid;
			string code;
			string length;
			string s;
			getline(ss, s, '.');
			getline(ss, lineid, '"');
			for (int i = 0; i < 4; ++i) {
				getline( ss, s, '=' );
			}
			getline( ss, code, ' ' );
			for (int i = 0; i < 6; ++i) {
				getline( ss, s, '=' );
			}
			getline( ss, length, '\n' );

			EDGES[hashTableEdge[lineid]].linecode = atoi(code.c_str());
			//NODE: Line data is bogus!!
			//EDGES[hashTableEdge[lineid]].length = atoi(length.c_str());

			//cout << lineid << " " << code << " " << length << endl;

		}
	}
	file3.close();

	ifstream file4(linecode_file.c_str());
	while (!file4.eof()) {
		getline(file4, line);
		if (line.size() > 0 && line[0] == 'N') {
			istringstream ss(line);
			string linecode;
			string numphases;
			string rvalues;
			string xvalues;
			vector<vector<double> > rmatrix;
			vector<vector<double> > xmatrix;
			string s;
			getline(ss, s, '.');
			getline(ss, linecode, '"');
			getline(ss, s, '=');
			getline(ss, numphases, ' ');
			getline( ss, s, '[' );
			getline( ss, rvalues, ']' );
			stringstream ssphases(numphases);
			int NumPhases;
			ssphases >> NumPhases;
			string row;
			istringstream ssrows(rvalues);
			for (int i = 0; i < NumPhases; ++i) { 
				rmatrix.push_back(vector<double>(NumPhases));
				getline(ssrows, row, '|');
				stringstream ssrow(row);
				for (int j = 0; j < NumPhases; ++j) {
					ssrow >> rmatrix[i][j];
					// CONVERT TO PER 1000 FEET
					//rmatrix[i][j] *= 5.28;
				}
			}
			getline( ss, s, '[' );
			getline( ss, xvalues, ']' );
			string row2;
			istringstream ssrows2(xvalues);
			for (int i = 0; i < NumPhases; ++i) { 
				xmatrix.push_back(vector<double>(NumPhases));
				getline(ssrows2, row2, '|');
				stringstream ssrow2(row2);
				for (int j = 0; j < NumPhases; ++j) {
					ssrow2 >> xmatrix[i][j];
					// CONVERT TO PER 1000 FEET
					//xmatrix[i][j] *= 5.28;
				}
			}
			pair<map<int, lineCodeData>::iterator, bool> ptr = LINECODES.insert(
				pair<int, lineCodeData>(atoi(linecode.c_str()),
				lineCodeData(atoi(linecode.c_str()), NumPhases, rmatrix, xmatrix)));



			/*for (int i = 0; i < NumPhases; ++i) {
				for (int j = 0; j < NumPhases; ++j) {
					cout << xmatrix[i][j] << " ";
				}
			}
			cout << endl;*/
		}
	}
	file4.close();

	/*ifstream file2(transformer_file.c_str());
	while (!file2.eof()) {
		getline(file2, line);
		if (line.size() > 0) {
			istringstream ss(line);
			string s;
			getline( ss, s, '.' );
			getline( ss, s, '"' );
			EDGES[hashTableEdge[s]].istransformer = true;
			EDGES[hashTableEdge[s]].linecode = -1;
		}
	}
	file2.close();

	// TRANSFORMERS
	vector<vector<double> > rtransformer;
	vector<vector<double> > xtransformer;

	for (int i = 0; i < 3; ++i) { 
		rtransformer.push_back(vector<double>(3, 0.0));
		xtransformer.push_back(vector<double>(3, 0.0));
	}
	pair<map<int, lineCodeData>::iterator, bool> ptr = LINECODES.insert(
		pair<int, lineCodeData>(-1, lineCodeData(-1, 3, rtransformer, xtransformer)));*/

	/*for (map<int, lineCodeData>::iterator it = LINECODES.begin(); it != LINECODES.end(); ++it) {
		cout << "LineCode " << it->first << " NumPhases " << it->second.numphases << endl;
	}*/

}	
