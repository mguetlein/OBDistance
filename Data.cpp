/*
 * Data.cpp
 *
 *  Created on: 28.07.2009
 *      Author: martin
 */

#include "Data.h"

using namespace std;
using namespace OpenBabel;

/*
Data::Data(char* smi_file, char* fragment_file, bool aromaticity) {

	this->aromaticity = aromaticity;

	readSmiles(smi_file);
	readSmarts(fragment_file);
}*/

Data::Data(char* smi_file, char* fragment_file, bool aromaticity, bool enable_3d, char* activity_file, char* dist_file) {

	this->aromaticity = aromaticity;
	this->enable_3d = enable_3d;

	readSmiles(smi_file);
	readSmarts(fragment_file);
	if (activity_file)
		readActivity(activity_file);
	if (dist_file)
		readDist(dist_file);
}

Data::~Data() {
	// TODO Auto-generated destructor stub
}


void Data::readSmiles(char * smiles_file) {
	string line;
	string tmp_field;
	string id;
	string smi;

	//vector<string> ids;
	//vector<string> smiles;


	//	string inchi;
	//	vector<string> ids, dup_ids;
	//	vector<string>::iterator dup_id;
	//	vector<string> inchis;
	//	vector<string>::iterator dup_inchi;

	//MolRef mol_ptr;
	int line_nr = 0;

	ifstream input;
	input.open(smiles_file);

	if (!input) {
		cerr << "Cannot open " << smiles_file << endl;
		//out->print_err();
		exit(1);
	}

	cerr << "Reading smiles from " << smiles_file << endl;
	//out->print_err();

	line_nr = 0;

	while (getline(input, line)) {

		istringstream iss(line);

		int field_nr = 0;
		id = -1;
		smi = "";

		while (getline(iss, tmp_field, '\t')) { // split at tabs

			if (field_nr == 0)
				id = tmp_field;
			else if (field_nr == 1)
				smi = tmp_field;

			field_nr++;
		}

		if (smi.length() == 0)
			cerr << "WARNING: empty smiles, line-nr: "<<line_nr<<", line: '"<<line<<"'"<<endl;

		ids.push_back(id);
		smiles.push_back(smi);

		// ID
		//dup_id = find(ids.begin(), ids.end(), id);

		/*if (dup_id == ids.end())
		 ids.push_back(id);

		 else {
		 *out << id << " (line " << line_nr << ") is not a unique ID ... exiting.\n";
		 out->print_err();
		 exit(1);
		 }*/

		// SMILES
		//remove_dos_cr(&smi);

		//mol_ptr = new FeatMol<MolType,FeatureType,ActivityType>(line_nr, id, smi,out);

		//inchi = mol_ptr->get_inchi();

		// do not print duplicate warning for empty smiles
		// (warning already printed in FeatMol constructor)
		/*if (inchi.size()>0)
		 {
		 dup_inchi = find(inchis.begin(), inchis.end(), inchi);

		 if (dup_inchi == inchis.end())
		 inchis.push_back(inchi);
		 else {
		 dup_ids =  this->get_idfrominchi(inchi);
		 *out << "Compounds " << id ;
		 for (dup_id=dup_ids.begin();dup_id!=dup_ids.end();dup_id++) {
		 *out << " and " << *dup_id;
		 }
		 *out << " have identical structures.\n";
		 out->print_err();
		 }
		 }

		 compounds.push_back(mol_ptr);
		 */
		line_nr++;
	}

	cerr << num_smiles() << " smiles"<< endl;

	input.close();
}

void Data::readSmarts(char * fragment_file) {
	string line;
	string tmp_field;
	string smart;
	vector<int> indices;

	int line_nr = 0;

	ifstream input;
	input.open(fragment_file);

	if (!input) {
		cerr << "Cannot open " << fragment_file << endl;
		//out->print_err();
		exit(1);
	}

	cerr << "Reading smarts from " << fragment_file << endl;
	//out->print_err();

	line_nr = 0;

	while (getline(input, line)) {

		istringstream iss(line);

		int field_nr = 0;
		smart = -1;
		indices.clear();

		while (getline(iss, tmp_field, '\t') && field_nr < 2) { // split at tabs

			if (field_nr==0)
				smart = tmp_field;
			else {

				size_t startpos = 0, endpos = 0;
				vector<string> tokens;

				for (;;) {

					startpos = tmp_field.find_first_not_of(" \t\n",startpos);
					endpos   = tmp_field.find_first_of(" \t\n",startpos);

					if (endpos < tmp_field.size() && startpos <= tmp_field.size())
						tokens.push_back(tmp_field.substr(startpos,endpos-startpos));
					else
						break;

					startpos = endpos + 1;

				}

				for (unsigned int i = 0 ; i < tokens.size() ; i++ ) {

					if ( tokens[i] == "[" )
						continue;
					else if ( tokens[i] == "]" )
						break;

					int comp_nr = atoi(tokens[i].c_str());

					indices.push_back(comp_nr);
					//feat_ptr->add_match(comp_nr);	// comp_nr is line number
					//this->get_compound(comp_nr)->add_feature(feat_ptr);

				}
			}
			field_nr++;
		}

		smarts.push_back(smart);
		smarts_occurences.push_back(indices);
		line_nr++;
	}

	cerr << num_smarts() << " smarts"<< endl;

	input.close();
}

void Data::readActivity(char * activity_file) {
	string line;
	string tmp_field;
	string id;
	string act;

	int line_nr = 0;

	ifstream input;
	input.open(activity_file);

	if (!input) {
		cerr << "Cannot open " << activity_file << endl;
		//out->print_err();
		exit(1);
	}

	cerr << "Reading activity from " << activity_file << endl;
	//out->print_err();

	line_nr = 0;

	num_actives = 0;
	num_inactives = 0;

	while (getline(input, line)) {

		istringstream iss(line);

		int field_nr = 0;
		id = -1;
		act = -1;

		while (getline(iss, tmp_field, '\t')) { // split at tabs

			if (field_nr == 0)
				id = tmp_field;
			else if (field_nr == 2)
				act = tmp_field;

			field_nr++;
		}

		if (act == "1")
			num_actives++;
		else
			num_inactives++;
		activity.push_back(act == "1");
		line_nr++;
	}

	cerr << num_actives << "/" << num_inactives << " active/inactive structures (smiles)" << endl;
	cerr << (num_actives / (double) num_inactives) << " active smiles ratio" << endl;

	input.close();
}


void Data::readDist(char * dist_file) {
	string line;
	string tmp_field;
	int idx1;
	int idx2;

	int line_nr = 0;

	ifstream input;
	input.open(dist_file);

	if (!input) {
		cerr << "Cannot open " << dist_file << endl;
		//out->print_err();
		exit(1);
	}

	cerr << "Reading distance pairs from " << dist_file << endl;
	//out->print_err();

	line_nr = 0;

	while (getline(input, line)) {

		istringstream iss(line);

		idx1 = -1;
		idx2 = -1;

		while (getline(iss, tmp_field, ']')) { // split at tabs

			tmp_field = tmp_field.substr(1,tmp_field.length()-1);
			int index = tmp_field.find_first_of(';',0);
//			cerr << tmp_field.substr(0,index).c_str() << endl;
			idx1 = atoi(tmp_field.substr(0,index).c_str());
//			cerr << tmp_field.substr(index+1,tmp_field.length()-(index+1)).c_str() << endl;
			idx2 = atoi(tmp_field.substr(index+1,tmp_field.length()-(index+1)).c_str());
			 break;
		}

//		cerr << idx1 << " " << idx2 << endl;
		int * indices = new int[2];
		indices[0] = idx1;
		indices[1] = idx2;
		dist.push_back(indices);
		line_nr++;
	}

	cerr << num_distance_pairs() << " distance pairs"<< endl;

	input.close();
}



OBMol * Data::get_mol(int smiles_index)
{
	if ( mols.find(smiles_index) == mols.end())
	{
		OBMol * mol = new OBMol();
		OBConversion obconversion;
		obconversion.SetInFormat("smiles");
		obconversion.ReadString(mol, smiles.at(smiles_index));

		//cerr << "unsetting aromatic: "<< smiles.at(smiles_index) << "\n";
		if (!aromaticity)
		{
			OBAtom a;
			FOR_ATOMS_OF_MOL(a,*mol)
				a->UnsetAromatic();
			OBBond b;
			FOR_ATOMS_OF_MOL(b,*mol)
				b->UnsetAromatic();
			mol->SetAromaticPerceived();
		}

		mols[smiles_index] = mol;
	}

	return mols[smiles_index];
}

vector<vector <int> > * Data::get_matches(int smiles_index, int smarts_index)
{
	return get_matches(smiles_index, smarts_index, false);
}

vector<vector <int> > * Data::get_matches(int smiles_index, int smarts_index, bool no_lookup)
{
	char buffer[255];
	sprintf(buffer, "%d_%d",smiles_index,smarts_index);
	string * key = new string(buffer);

	if ( matches.find(*key) == matches.end() || no_lookup )
	{
		OBMol * mol = get_mol( smiles_index );
		OBSmartsPattern smartsPattern;
		smartsPattern.Init( *get_smarts( smarts_index ) );
		smartsPattern.Match( *mol );
		vector<vector <int> > map = smartsPattern.GetMapList();

		if (no_lookup)
		{
			cerr << "smiles "<< smiles.at(smiles_index) << "\n";
			cerr << "formula "<< mol->GetFormula() << "\n";
			cerr << "smarts "<< *get_smarts( smarts_index ) << "\n";
			cerr << "matches " << map.size()<< "\n";
		}

		vector<vector <int> > * disjunced_map = new vector<vector <int> >();
		for (unsigned int i = 0 ; i < map.size(); i++ )
		{
//			if (no_lookup)
//				cerr << "check for match "<< i << "\n";
//
//			//check for match
//			bool match = false;
//			for (unsigned int j = 0 ; j < disjunced_map->size(); j++ )
//			{
//				for (unsigned int k = 0 ; k < map[i].size(); k++ )
//				{
//					for (unsigned int l = 0 ; l < disjunced_map->at(j).size(); l++ )
//					{
//						if (map[i][k] == disjunced_map->at(j)[l])
//						{
//							//cerr << "match - ommit "<< i << "\n";
//							match = true;
//							break;
//						}
//					}
//					if (match)
//						break;
//				}
//				if(match)
//					break;
//			}
//
//			if (!match)
				disjunced_map->push_back(map[i]);
		}
		//cerr << "'" << *key << "'" << endl;

		matches[*key] = disjunced_map;
	}

	return matches[*key];

}
