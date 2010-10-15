/*
 * Data.h
 *
 *  Created on: 28.07.2009
 *      Author: martin
 */

#ifndef DATA_H_
#include "OBDistance.h"
#define DATA_H_

using namespace std;
using namespace OpenBabel;

class Data {

public:

	Data(char* smi_file, char* fragment_file, bool aromaticity, bool enable_3d, char* activity_file, char* dist_file);
	virtual ~Data();

	void readSmiles(char * smiles_file);
	void readSmarts(char * fragment_file);
	void readActivity(char * activity_file);
	void readDist(char * dist_file);

	string * get_id(int smiles_index)
	{
		return &ids[smiles_index];
	};

	unsigned int num_smiles()
	{
		return smiles.size();
	}
	string * get_smiles(int smiles_index)
	{
		return &smiles[smiles_index];
	};

	bool is_active(int smiles_index)
	{
		return activity[smiles_index];
	};

	unsigned int num_smarts()
	{
		return smarts.size();
	}
	string * get_smarts(int smarts_index)
	{
		return &smarts[smarts_index];
	};

	vector<int> * get_smarts_occurences(int smarts_index)
	{
		return &smarts_occurences[smarts_index];
	}

	OBMol * get_mol(int smiles_index);
	vector<vector <int> > * get_matches(int smiles_index, int smarts_index );
	vector<vector <int> > * get_matches(int smiles_index, int smarts_index, bool no_lookup);


	int * get_distance_pair(int dist_index)
	{
		return dist[dist_index];
	}
	unsigned int num_distance_pairs()
	{
		return dist.size();
	}

	bool is_3d_enabled()
	{
		return enable_3d;
	}

	unsigned int get_num_actives()
	{
		return num_actives;
	}

	unsigned int get_num_inactives()
	{
		return num_inactives;
	}

private:

	bool aromaticity;
	bool enable_3d;

	vector<string> ids;
	vector<string> smiles;
	vector<bool> activity;
	vector<string> smarts;
	vector<vector <int> > smarts_occurences;
	vector<int *> dist;
	unsigned int num_actives;
	unsigned int num_inactives;

	map<int, OBMol *> mols;
	map<string, vector<vector <int> > *> matches;
};

#endif /* DATA_H_ */
