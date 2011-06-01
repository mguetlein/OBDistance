/*
 * CheckFragments.cpp
 *
 *  Created on: 27.07.2009
 *      Author: martin
 */

#include "CheckSmarts.h"

using namespace std;
using namespace OpenBabel;


CheckSmarts::CheckSmarts(Data * data) {

	this->data = data;
}

CheckSmarts::~CheckSmarts() {
}



void CheckSmarts::checkSmarts() {

	OBSmartsPattern smartsPattern;

	for (unsigned int i = 0; i < data->num_smarts(); i++) {

		cout << *data->get_smarts(i) << "\t[ ";

		for (unsigned int j = 0; j < data->num_smiles(); j++) {

			smartsPattern.Init(*data->get_smarts(i));

			if (smartsPattern.Match(*data->get_mol(j), true))
				cout << data->get_id(j) << " ";
		}

		cout << "]\n";
	}
}

void CheckSmarts::validateSmarts() {

	OBSmartsPattern smartsPattern;

	int matches = 0;
	int noMatches = 0;
	unsigned int min_freq = 1000;

	for (unsigned int j = 0; j < data->num_smarts(); j++) {

		vector<int > * occ = data->get_smarts_occurences(j);

		if (occ->size() < min_freq)
			min_freq = occ->size();

		for (unsigned int i = 0; i < occ->size(); i++)
		{
			vector<vector <int> > * map1 = data->get_matches_from_id(occ->at(i), j);
			if (map1->size() < 1)
			{
				cerr << "no match, smarts: " << *data->get_smarts(j)<< ",  smiles: "<< *data->get_smiles_from_id(occ->at(i)) << endl;
				noMatches++;
			}
			else
				matches++;
		}
	}

	cout << "validated matches: "<<matches<<endl;
	cout << "no matches: "<<noMatches<<endl;
	cout << "minimum frequency: "<<min_freq<<endl;
}
