#include <iostream>
#include "rapidxml-1.13/rapidxml.hpp"
#include "rapidxml-1.13/rapidxml_print.hpp"
//#include "rapidxml-1.13/rapidxml_iterators.hpp"
//#include "rapidxml-1.13/rapidxml_utils.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string.h>
#include <stdlib.h>

using namespace rapidxml;
using namespace std;

void AlgorithmData(xml_node<> *pRoot, double &LoadMet, double &CriticalLoadMet, double &PhaseVariation) {
	for (xml_node<> *pNode = pRoot->first_node(); pNode; pNode = pNode->next_sibling()) {
		string str = pNode->name();
		if (str.compare("LoadMet") == 0) {
			LoadMet = atof(pNode->value());
		} else if (str.compare("CriticalLoadMet") == 0) {
			CriticalLoadMet = atof(pNode->value());
		} else if (str.compare("PhaseVariation") == 0) {
			PhaseVariation = atof(pNode->value());
		}
	}
}

void readModification (xml_node<> *pNode, xml_node<> *pState, string &name, string &value) {
	xml_node<> *pNameNode = pNode->first_node("Identifiers");
	pNameNode = pNameNode->first_node("Identifier");
	pNameNode = pNameNode->first_node("Value");
	pState = pState->next_sibling();
	name = pNameNode->value();
	value = pState->value();
	//DISABLED_name.push_back(name);
	//name = pState->value();
	//DISABLED_data.push_back(name.compare("true") == 0);
}

template<class T> struct dataNode {
public:
	dataNode() {}
	dataNode(string name, T x) {id = name; data = x;}
	//string type;
	string id;
	T data;
};

void SlaveModificationsData(xml_node<> *pRoot, 
	vector<dataNode<bool> > &DISABLED,
	vector<dataNode<bool> > &HARDENED_DISABLED) {
	for (xml_node<> *pNode = pRoot->first_node(); pNode; pNode = pNode->next_sibling()) {
		xml_node<> *pSubNode = pNode->first_node("Attribute");
		xml_node<> *pState = pSubNode->first_node("Key");
		string str = pState->value();
		if (str.compare("DISABLED") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			// MAGIC LINES
			if (false) { //name[0] == 'n' || name[0] == 'o') {
				DISABLED.push_back(dataNode<bool>(name,false)); 
			} else {
				DISABLED.push_back(dataNode<bool>(name,value.compare("true") == 0)); 
			}
		} else if (str.compare("HARDENED_DISABLED") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			// MAGIC LINES
			if (false) { //name[0] == 'n' || name[0] == 'o') {
				HARDENED_DISABLED.push_back(dataNode<bool>(name,false)); 
			} else {
				HARDENED_DISABLED.push_back(dataNode<bool>(name,value.compare("true") == 0)); 
			}
		} else {
			cout << "WHOAAA!!!" << endl;
		}
	}
}

void MasterModificationsData(xml_node<> *pRoot, 
	vector<dataNode<double> > &LINE_CONSTRUCTION_COST,
	vector<dataNode<double> > &MICROGRID_COST,
	vector<dataNode<double> > &MICROGRID_FIXED_COST,
	vector<dataNode<bool> > &IS_CRITICAL_LOAD,
	vector<dataNode<double> > &MAX_MICROGRID,
	vector<dataNode<double> > &HARDEN_COST,
	vector<dataNode<double> > &LINE_SWITCH_COST,
	vector<dataNode<string> > &BUILD_SET,
	vector<dataNode<string> > &HARDEN_SET) {
	for (xml_node<> *pNode = pRoot->first_node(); pNode; pNode = pNode->next_sibling()) {
		xml_node<> *pSubNode = pNode->first_node("Attribute");
		xml_node<> *pState = pSubNode->first_node("Key");
		string str = pState->value();
		if (str.compare("LINE_CONSTRUCTION_COST") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			// MAGIC LINES
			//LINE_CONSTRUCTION_COST.push_back(dataNode<double>(name, 0.0));
			LINE_CONSTRUCTION_COST.push_back(dataNode<double>(name, atof(value.c_str())));
		} else if (str.compare("MICROGRID_COST") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			MICROGRID_COST.push_back(dataNode<double>(name,atof(value.c_str())));
		} else if (str.compare("MICROGRID_FIXED_COST") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			MICROGRID_FIXED_COST.push_back(dataNode<double>(name,atof(value.c_str())));
		} else if (str.compare("IS_CRITICAL_LOAD") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			IS_CRITICAL_LOAD.push_back(dataNode<bool>(name,value.compare("true") == 0)); 
		} else if (str.compare("MAX_MICROGRID") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			MAX_MICROGRID.push_back(dataNode<double>(name,atof(value.c_str())));
		} else if (str.compare("HARDEN_COST") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			HARDEN_COST.push_back(dataNode<double>(name,atof(value.c_str())));
		} else if (str.compare("LINE_SWITCH_COST") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			LINE_SWITCH_COST.push_back(dataNode<double>(name,atof(value.c_str())));
		} else if (str.compare("BUILD_SET") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			BUILD_SET.push_back(dataNode<string>(name,value));
		} else if (str.compare("HARDEN_SET") == 0) {
			string name;
			string value;
			readModification(pNode, pState, name, value);
			HARDEN_SET.push_back(dataNode<string>(name,value));
		} else {
			cout << "WHOAAA!!!" << endl;
		}
	}
}

void readSLAVEData(string input_file, 
	vector<dataNode<bool> > &DISABLED,
	vector<dataNode<bool> > &HARDENED_DISABLED) {

	xml_document<> doc;
	ifstream file(input_file.c_str());

	stringstream buffer;
	buffer << file.rdbuf();
	file.close();
	string content(buffer.str());
	doc.parse<0>(&content[0]);

	xml_node<> *pRoot = doc.first_node();

	for (xml_node<> *pNode = pRoot->first_node(); pNode; pNode = pNode->next_sibling()) {
		string str = pNode->name();
		if (str.compare("Modifications") == 0) {
			SlaveModificationsData(pNode, DISABLED, HARDENED_DISABLED);
		}
	}

	/*for (size_t i = 0; i < DISABLED.size(); ++i) {
		cout << DISABLED[i].id << " " << DISABLED[i].data << endl;
	}*/
}

// SAMPLE USE:
void readMASTERData(string input_file, 
	vector<dataNode<double> > &LINE_CONSTRUCTION_COST,
	vector<dataNode<double> > &MICROGRID_COST,
	vector<dataNode<double> > &MICROGRID_FIXED_COST,
	vector<dataNode<bool> > &IS_CRITICAL_LOAD,
	vector<dataNode<double> > &MAX_MICROGRID,
	vector<dataNode<double> > &HARDEN_COST,
	vector<dataNode<double> > &LINE_SWITCH_COST,
	vector<dataNode<string> > &BUILD_SET,
	vector<dataNode<string> > &HARDEN_SET,
	double &LoadMet,
	double &CriticalLoadMet,
	double &PhaseVariation) {

	xml_document<> doc;
	ifstream file(input_file.c_str());

	stringstream buffer;
	buffer << file.rdbuf();
	file.close();
	string content(buffer.str());
	doc.parse<0>(&content[0]);

	xml_node<> *pRoot = doc.first_node();

	for (xml_node<> *pNode = pRoot->first_node(); pNode; pNode = pNode->next_sibling()) {
		string str = pNode->name();
		if (str.compare("Algorithm") == 0) {
			AlgorithmData(pNode, LoadMet, CriticalLoadMet, PhaseVariation);
			//cout << LoadMet << endl;
			//cout << CriticalLoadMet << endl;
			//cout << PhaseVariation << endl;
		} else if (str.compare("Modifications") == 0) {
			MasterModificationsData(pNode, LINE_CONSTRUCTION_COST, MICROGRID_COST,
					MICROGRID_FIXED_COST, IS_CRITICAL_LOAD, MAX_MICROGRID, HARDEN_COST, LINE_SWITCH_COST, BUILD_SET, HARDEN_SET);
		}
	}

	/*for (size_t i = 0; i < DISABLED.size(); ++i) {
		cout << DISABLED[i].id << " " << DISABLED[i].data << endl;
	}*/
}
