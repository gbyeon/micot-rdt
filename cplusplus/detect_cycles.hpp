// AUTHOR: EMRE YAMANGIL, 2014
// Contact: emreyamangil@gmail.com, emreyamangil@lanl.gov
// Description:
// 	Following code performs a biconnected component decomposition of an undirected graph and uses a DFS variant to find all cycles within each biconnected component. 

#include <iostream>
#include <vector>
#include <list>
#include <stack>
#include <fstream>
#include <set>
#include <algorithm> //for std::rotate(myvector.begin(),myvector.begin()+i,myvector.end()); // i is the rotate index (i = 1 will rotate once to right)
using namespace std;


// Graph vertex structure
class graph_vertex {
public:
	graph_vertex() {visited = false; included = false; parent = -1;};
	int label;
	string id;
	list<int> AdjList;
	list<string> EdgeID;
	bool visited;
	bool included;
	int parent;
	int depth;
	int low;
	int pred;	
};

// Data is written in:
// N M (N number of vertices, M number of edges)
// vi1 vi2 (where edge ei = (vi1 vi2)
/*void read_data (vector<graph_vertex> &G) {
	ifstream file("graph.txt");
	int N, M;
	int u, v;
	file >> N >> M;
	G.resize(N);
	for (int i = 0; i < M; ++i) {
		file >> u >> v;
		G[u].AdjList.push_back(v);
		G[v].AdjList.push_back(u);
	}
}*/

void OutputComponent (stack<vector<int> > &S, vector<int> &e, vector<list<vector<int> > > &C) {
	list<vector<int> > Comp;
	vector<int> f;
	do {
		f = S.top();
		//cout << f[0] << f[1] << endl;
		Comp.push_back(f);
		S.pop();
	} while (f[0] != e[0] || f[1] != e[1]);
	C.push_back(Comp);
}

int min (int a, int b) {
	return a < b ? a : b;
}

void DFS_visit(vector<graph_vertex> &G, stack<vector<int> > &S, int u, int &depth, vector<list<vector<int> > > &C) {
	depth = depth+1;
	G[u].visited = true;
	G[u].depth = depth;
	G[u].low = depth;
	for (list<int>::iterator it = G[u].AdjList.begin(); it != G[u].AdjList.end(); ++it) {
		int v = *it;
		if (!(G[v].visited)) {
			//cout << v << endl;
			vector<int> e(2);
			e[0] = u; e[1] = v;
			S.push(e);
			G[v].parent = u;
			DFS_visit(G, S, v, depth, C);
			if (G[v].low >= G[u].depth) {
				OutputComponent(S, e, C);
			}
			G[u].low = min(G[u].low, G[v].low);
		} else if ((G[u].parent != v) && (G[v].depth < G[u].depth)) {
			// uv is a back edge from u to its ancestor v
			vector<int> e(2);
			e[0] = u; e[1] = v;
			S.push(e);
			G[u].low = min(G[u].low, G[v].depth);
		}
	}
}

void DFS(vector<graph_vertex> &G, stack<vector<int> > &S, int source, vector<list<vector<int> > > &C) {
	int depth = 0;
	int N = G.size();
	for (int v = 0; v < N; ++v) {
		if (!G[v].visited) {
			DFS_visit (G, S, v, depth, C); 
		}
	}
}

void subGraphDFS_visit(vector<graph_vertex> &G, list<vector<int> > &S, int u, int &depth, set<int> &V) {
	depth = depth+1;
	G[u].visited = true;
	G[u].depth = depth;
	G[u].low = depth;
	for (list<int>::iterator it = G[u].AdjList.begin(); it != G[u].AdjList.end(); ++it) {
		int v = *it;
		if (!(G[v].visited) && V.find(v) != V.end()) {
			G[v].parent = u;
			subGraphDFS_visit(G, S, v, depth, V);
		} else if ((G[u].parent != v) && V.find(v) != V.end()) {
			// uv is a back edge
			vector<int> e(2);
			e[0] = u; e[1] = v;
			S.push_back(e);
		}
	}
}

void subGraphDFS(vector<graph_vertex> &G, list<vector<int> > &S, set<int> &V) {
	int depth = 0;
	int N = G.size();
	for (int i = 0; i < N; ++i) {
		G[i].visited = false;
	}
	for (set<int>::iterator it = V.begin(); it != V.end(); ++it) {
		if (!G[*it].visited) {
			subGraphDFS_visit (G, S, *it, depth, V); 
		}
	}
}

void GenerateSubGraph(list<vector<int> > &C, set<int> &V) {
	for (list<vector<int> >::iterator it = C.begin(); it != C.end(); ++it) {
		V.insert((*it)[0]);
		V.insert((*it)[1]);
	}	
}

bool visited(vector<int> &path, int v) {
	for (int i = 0; i < path.size(); ++i) {
		if (path[i] == v) {
			return true;
		}
	}
	return false;	
}

int findMinIndex(vector<int> &path) {
	int min_val = 1e+5;
	int min_ind = -1;
	for (int i = 0; i < path.size(); ++i) {
		if (path[i] < min_val) {
			min_val = path[i];
			min_ind = i;	
		}
	}
}

bool equalPath (const vector<int> &p_1, const vector<int> &p_2) {
	int n = p_1.size();
	for (int j = 0; j < n; ++j) {
		if (p_1[j] != p_2[j]) {
			return false;
		}
	}
	return true;
}

bool isNewPath (vector<vector<int> > &P, vector<int> path) {
	int n = path.size();
	for (int i = 0; i < P.size(); ++i) {
		if (n == P[i].size()) {
			if (equalPath(path, P[i])) {
				return false;
			}
		}
	}
	return true;
}

void insertPath(vector<vector<int> > &P, vector<int> path) {
	//int min_ind = findMinIndex(path);
	std::sort(path.begin(), path.end());
	if (isNewPath(P, path)) {
		P.push_back(path);
	}
}

void findAllstPaths(vector<graph_vertex> &G, vector<vector<int> > &P, vector<int> path, int u, int t, set<int> &V) {
	for (list<int>::iterator it = G[u].AdjList.begin(); it != G[u].AdjList.end(); ++it) {
		int v = *it;
		if (v == t && path.size() == 1) {
			continue;
		}
		if (!visited(path, v) && V.find(v) != V.end()) {
			if (v == t) {
				path.push_back(v);
				//P.push_back(path);
				insertPath(P, path);
				path.pop_back();
			} else {
				path.push_back(v);
				findAllstPaths (G, P, path, v, t, V);
				path.pop_back();
			}
		}
	}
}

void ComponentVisit(vector<graph_vertex> &G, int i, vector<int> &newComp) {
	G[i].visited = true;
	for (list<int>::iterator it = G[i].AdjList.begin(); it != G[i].AdjList.end(); ++it) {
		if (!G[*it].visited) {
			newComp.push_back(*it);
			ComponentVisit(G, *it, newComp);
		}
	}
}

void findComponents(vector<graph_vertex> &G, list<vector<int> > &components) {
	// Initialize DFS
	for (int i = 0; i < G.size(); ++i) {
		G[i].visited = false;
	}
	// Compute components
	for (int i = 0; i < G.size(); ++i) {
		if (!G[i].visited) {
			vector<int> newComp;
			newComp.push_back(i);
			ComponentVisit(G, i, newComp);
			components.push_back(newComp);
		}
	}
}

void subComponentVisit(vector<graph_vertex> &G, int source, int i, list<vector<int> > &edges) {
	G[i].visited = true;
	for (list<int>::iterator it = G[i].AdjList.begin(); 
		it != G[i].AdjList.end(); ++it) {
		int v = *it;
		if (!G[v].visited) {
			G[v].pred = i;
			if (G[v].included) {
				int pred = i;
				int curr = v;
				while (curr != source) {
					vector<int> edge(2);
					edge[0] = pred;
					edge[1] = curr;
					std::sort(edge.begin(), edge.end());
					curr = pred;
					pred = G[pred].pred;
					edges.push_back(edge);
				}
			}
			subComponentVisit(G, source, v, edges);
		}
	}
}

void findSubComponent(vector<graph_vertex> &G, int source, list<vector<int> > &edges) {
	// Initialize DFS
	for (int i = 0; i < G.size(); ++i) {
		G[i].visited = false;
	}
	subComponentVisit(G, source, source, edges);
}

void detectCycles(vector<graph_vertex> &G, vector<vector<int> > &cycles) {
	//vector<graph_vertex> G;
	//read_data(G);
	//vector<vector<int> > cycles;

	// COMPONENTS
	vector<list<vector<int> > > C;
	stack<vector<int> > S;
	int source = 0;
	DFS(G, S, source, C);
	/*for (int i = 0; i < C.size(); ++i) {
		cout << "NEW COMPONENT!" << endl;
		for (list<vector<int> >::iterator it = C[i].begin(); it != C[i].end(); ++it) {
			cout << "(" << (*it)[0] << "," << (*it)[1] << ") ";
		}
		cout << endl;
	}*/
		
	for (int i = 0; i < C.size(); ++i) {
		set<int> V;
		GenerateSubGraph(C[i], V);
		list<vector<int> > E;
		subGraphDFS(G, E, V);
		/*for (set<int>::iterator it = V.begin(); it != V.end(); ++it) {
			cout << *it;
		}
		cout << endl;*/
		for (list<vector<int> >::iterator it = E.begin(); it != E.end(); ++it) {
			vector<int> path;
			//cout << "BACK EDGE " << (*it)[0] << (*it)[1] << endl;
			path.push_back((*it)[0]);
			findAllstPaths(G, cycles, path, (*it)[0], (*it)[1], V);
		}
	}
	for (int i = 0; i < cycles.size(); ++i) {
		cout << "CYCLE #" << i << ": "; 		
		for (int j = 0; j < cycles[i].size(); ++j) {
			cout << cycles[i][j] << ", ";
		}
		cout << endl;	
	}
	
}
