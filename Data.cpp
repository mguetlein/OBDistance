/*
 * Data.cpp
 *
 *  Created on: 28.07.2009
 *      Author: martin
 */

#include <boost/regex.hpp>

#include "Data.h"

using namespace std;
using namespace OpenBabel;

extern bool DEBUG_OUT;
extern bool CACHE_FRAGMENTS;
/*
Data::Data(char* smi_file, char* fragment_file, bool aromaticity) {

	this->aromaticity = aromaticity;

	readSmiles(smi_file);
	readSmarts(fragment_file);
}*/

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

string noHydrogenFormula( OBMol * mol )
{
//	 string formula = mol->GetFormula();
//     return tr1::regex_replace(formula, rx("H"), "");

     //std::locale::global(std::locale("German"));
	 boost::regex expr("H[0-9]*");
	 std::string fmt("\\2\\1");
	 std::string res = boost::regex_replace(mol->GetFormula(), expr, fmt);
	 return res;

}


Data::Data(char* smi_file, char* fragment_file, bool aromaticity, bool enable_3d,
		char* activity_file, char* dist_file, char* sdf_file) { //, char* sdf_smiles_file) {

	this->aromaticity = aromaticity;
	this->enable_3d = enable_3d;

	readSmiles(smi_file);
	readSmarts(fragment_file);
	if (activity_file)
		readActivity(activity_file);
	if (dist_file)
		readDist(dist_file);
	if (sdf_file)
		readMolsFromSdf(sdf_file); //, sdf_smiles_file);
}

Data::~Data() {
	// TODO Auto-generated destructor stub
}

void Data::readMolsFromSdf(char * sdf_file) { //, char * sdf_smiles_file) {

	cerr << "reading mols from sdf file: "<< sdf_file << endl; //", with according smiles file: " << sdf_smiles_file << endl;

//	vector<string> smiles;
//	string line;
//	ifstream myfile (sdf_smiles_file);
//	if (myfile.is_open())
//	{
//	  while ( myfile.good() )
//	  {
//	    getline (myfile,line);
//	    //cerr << line << endl;
//	    if (line.size()>0)
//	    	smiles.push_back(line);
//	  }
//	  myfile.close();
//	}
//    else
//    {
//    	cerr << "Unable to open file "<< sdf_smiles_file<<endl;
//    }

	OBConversion obconversion;
	obconversion.SetInFormat("sdf");
	obconversion.SetOutFormat("smiles");

//	OBConversion obconversion2;
//	obconversion2.SetInFormat("smiles");

	OBMol * mol = new OBMol();
	bool notatend = obconversion.ReadFile(mol,sdf_file);
	unsigned int count = 0;
	while (notatend)
	{


//		string inchi = obconversion.WriteString(mol,true);
//		if ( sdfMolsSmiles.find(inchi) != sdfMolsSmiles.end())
//		{
//			cerr << "duplicate inchi "+inchi << endl;
////			if ( abs(sdfMolsSmiles[inchi]->GetMolWt() - mol->GetMolWt() ) > 0.001 )
////			{
////				cerr << "diff mol weight "<<sdfMolsSmiles[inchi]->GetMolWt()<<" != "<<mol->GetMolWt()<<endl;
////				exit(1);
////			}
//		}
//		else
//		{
//			cerr << "<< "<<inchi<<endl;
//			if ( smiles.at(count) != obconversion.WriteString(mol,true) )
//			{
//				cerr << "wrong smiles at "<<count<<endl;

//				OBMol m;
//				obconversion2.ReadString(&m, smiles.at(count));
//
//				if (mol->GetFormula() != m.GetFormula())
//				{
//					cerr << "formula does not fit" << endl;
//					cerr << count;
//					cerr << smiles.at(count) << endl;
//					cerr << "formula from smiles " << m.GetFormula() << endl;
//					cerr << m.GetMolWt() << endl;
//					cerr << "formula from sdf    " << mol->GetFormula() << endl;
//					cerr << mol->GetMolWt() << endl;
//					exit(1);
//				}

				//cerr << count << endl;
				//cerr << smiles.at(count) << endl;
				//cerr << obconversion.WriteString(mol,true) << endl;
//				exit(1);
//			}
			//sdfMolsSmiles[smiles.at(count)] = mol;


//			cerr << mol->GetFormula() << endl;
//			sdfMolsSmiles[mol->GetFormula()] = mol;

			// for some pdb files the compound name is added to the smiles
			string smi = obconversion.WriteString(mol,true);
			vector<string> v;
			split(smi,'\t',v);
			if (v.size() > 0)
				smi = v[0];
			string formula = noHydrogenFormula(mol);

			int match = -1;
			for (unsigned int i = 0 ; i < smiles.size() ; i++ )
			{
				if (smi == smiles[i])
				{
					match = i;
					break;
				}
			}
			if (match == -1)
			{
				for (unsigned int i = 0 ; i < smiles.size() ; i++ )
				{
					OBMol m;
					OBConversion obconversion;
					obconversion.SetInFormat("smiles");
					obconversion.ReadString(&m, smiles.at(i));
					string f = noHydrogenFormula(&m);
					cerr << "    " << f << endl;
					if (f == formula)
					{
						if (match != -1)
						{
							cerr << "two formulas match" << endl;
							exit(1);
						}
						cerr << "    match" << endl;
						match = i;
					}
				}
			}
			if (match == -1)
				cerr << "no match in smiles file found" << endl;
			if (DEBUG_OUT || match==-1)
			{
				cerr << sdfMolsSmiles.size() << endl;
				cerr << "  " << smi << endl;
				cerr << "  " << formula << endl; // mol->GetFormula() << endl;
				cerr << "  " << mol->Has3D() << endl;

			}
//			if (match == -1)
//				exit(1);

			if (!mol->Has3D())
				cerr << "  warning: no 3D info for " << smi << endl;

			if (match != -1)
				sdfMolsSmiles[smiles[match]] = mol;

//			sdfMolsSmiles[smi] = mol;
			//sdfMolsFormula[noHydrogenFormula(mol)] = mol;
//		}
		mol = new OBMol();
		notatend = obconversion.Read(mol);
		count++;
	}
//	if (smiles.size() != count)
//	{
//		cerr << "size does matter: " << smiles.size() << " != " << count << endl;
//		exit(1);
//	}
	//exit(1);
	//delete mol;
}



void Data::readSmiles(char * smiles_file) {
	string line;
	string tmp_field;
	int id;
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
				id = atoi(tmp_field.c_str());
			else if (field_nr == 1)
				smi = tmp_field;

			field_nr++;
		}

		if (smi.length() == 0)
			cerr << "WARNING: empty smiles, line-nr: "<<line_nr<<", line: '"<<line<<"'"<<endl;

		ids.push_back(id);
		id_to_smiles_index[id] = smiles.size();
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
	int id;
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
			{
				id = atoi(tmp_field.c_str());
				if (id != ids[line_nr])
				{
					cerr << "order smiles & classes is not equal, implement id instead of line usage for class values" << endl;
					exit(1);
				}
			}
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

	if (activity.size() != smiles.size())
	{
		cerr << "Number of activties (" << activity.size() <<
			") does not match number of smiles (" << smiles.size() << ")" << endl;
		exit(1);
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

		if ( (min(idx1,idx2) < 0) || (max(idx1,idx2) >= (signed int)smarts.size()))
		{
			cerr << "Invalid smarts indices: " << tmp_field << endl;
			exit(1);
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
	if (DEBUG_OUT)
		cerr << "      get mol for smiles_index: " << smiles_index << ", smiles: " << smiles.at(smiles_index) << endl;

	if ( mols.find(smiles_index) == mols.end())
	{
		if (DEBUG_OUT)
			cerr << "        not yet created" << endl;

		OBMol m;
		OBConversion obconversion;
		obconversion.SetInFormat("smiles");
		obconversion.ReadString(&m, smiles.at(smiles_index));
		string formula = noHydrogenFormula(&m); // m.GetFormula();

		if (DEBUG_OUT)
			cerr << "        formula is " << formula << endl;

		OBMol * mol;

		if (sdfMolsSmiles.size() == 0)
		{
			mol = new OBMol();
			OBConversion obconversion;
			obconversion.SetInFormat("smiles");
			obconversion.ReadString(mol, smiles.at(smiles_index));

			if (this->enable_3d)
			{
				// build 3d coordinates
				if (DEBUG_OUT)
					cerr << "    build 3D ...";
				OBBuilder builder;
				builder.Build(*mol);
				if (DEBUG_OUT)
					cerr << "    done" << endl;
			}
		}
		else if (sdfMolsSmiles.size() > 0)
		{
			if (DEBUG_OUT)
				cerr << "        load from sdf"<< endl;

//			obconversion.SetOutFormat("inchi");
//			string inchi = obconversion.WriteString(&m,true);

			if ( sdfMolsSmiles.find(smiles.at(smiles_index)) == sdfMolsSmiles.end())
			//if ( sdfMolsSmiles.find(formula) == sdfMolsSmiles.end())
			{
//				bool formulaFound = false;
//				if ( sdfMolsFormula.find(formula) == sdfMolsFormula.end())
//				{
//					formulaFound = true;
//				}


				cerr << "could not find mol in sdf file" << endl;
				cerr << smiles_index << endl;
				cerr << smiles.at(smiles_index) << endl;
				cerr << formula << endl;
//				cerr << formulaFound << endl;
				exit(1);

//				mol = new OBMol();
//				OBConversion obconversion;
//				obconversion.SetInFormat("smiles");
//				obconversion.ReadString(mol, smiles.at(smiles_index));
			}
			else
			{
				mol = sdfMolsSmiles[smiles.at(smiles_index)];
				string form = noHydrogenFormula(mol);
				//mol = sdfMolsSmiles[formula];

				if (DEBUG_OUT)
				{
					cerr << "        found in sdf with equal smiles, checking formula"<< endl;
					//cerr << "        " << sdfMolsSmiles[smiles.at(smiles_index)]->GetFormula() << endl;
					//cerr << "        " << mol->GetFormula() << endl;
				}
				if (form != formula)
				{
					cerr << "formula does not fit" << endl;
					cerr << smiles_index << endl;
					cerr << smiles.at(smiles_index) << endl;
					cerr << formula << endl;
					cerr << form << endl;
					exit(1);
				}

				if (DEBUG_OUT)
				{
					cerr << "        formula is equal"<< endl;
					cerr << "        3d available "<< mol->Has3D() << endl;
				}
			}
		}

		//cerr << "unsetting aromatic: "<< smiles.at(smiles_index) << "\n";
		if (!aromaticity)
		{
			if (DEBUG_OUT)
				cerr << "        unset aromaticity"<< endl;

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

	if (DEBUG_OUT)
	{
		cerr << "      mol loading done"<< endl;
		cerr << "      3d available "<< mols[smiles_index]->Has3D() << endl;
	}
	return mols[smiles_index];
}

void Data::free_memory(int smiles_index, int smarts_index)
{
	unsigned long key = smiles_index * 10000000 + smarts_index;
	if ( matches.find(key) != matches.end())
	{
		if(DEBUG_OUT)
		{
//			map<unsigned long, vector<vector <int> > *>::iterator it;
	//		for ( it=matches.begin() ; it != matches.end(); it++ )
	//		  cerr << (*it).first << " => ?,  " << endl;
			cerr << "  free matches with key: " << key << "\n" << flush;
			cerr << "  num matches stored (before free): " << matches.size()<< "\n" << flush;
		}

		vector<vector <int> > * disjunced_map = matches[key];
		delete(disjunced_map);
		matches.erase(key);
	}
}

vector<vector <int> > * Data::get_matches(int smiles_index, int smarts_index)
{
	return get_matches(smiles_index, smarts_index, false);
}

vector<vector <int> > * Data::get_matches(int smiles_index, int smarts_index, bool no_lookup)
{
	unsigned long key = smiles_index * 10000000 + smarts_index;
	// example
	// smiles: 132, smarts: 13002
	// -> key: 1320013002

	if ( matches.find(key) == matches.end() || no_lookup )
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

		matches[key] = disjunced_map;
	}
	else if(!CACHE_FRAGMENTS && matches.size()>1)
	{
		cerr << "WTF" << endl;
		exit(1);
	}

	return matches[key];

}
