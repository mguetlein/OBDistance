/*
 * MineDistances.cpp
 *
 *  Created on: 28.07.2009
 *      Author: martin
 */

#include "MineDistances.h"

using namespace std;
using namespace OpenBabel;

extern bool DEBUG_OUT;
extern bool CACHE_FRAGMENTS;

MineDistances::MineDistances(Data * data,
		double rel_min_freq_per_class, double freq_ratio_tolerance) {
		//unsigned int min_frequency, unsigned int min_frequency_per_class) {

	this->data = data;

	if ( rel_min_freq_per_class > 0 )
	{
		abs_freq_active = max(1,(int)(data->get_num_actives() * rel_min_freq_per_class));
		abs_freq_inactive = max(1,(int)(data->get_num_inactives() * rel_min_freq_per_class));
		cerr << abs_freq_active << "/" << abs_freq_inactive << " minimum frequency in active/inactive structures" << endl;
		double active_ratio = data->get_num_actives() / (double) data->get_num_inactives();
		freq_ratio_min = active_ratio - active_ratio * freq_ratio_tolerance;
		freq_ratio_max = active_ratio + active_ratio * freq_ratio_tolerance;
		cerr << "[" << freq_ratio_min << " - " << freq_ratio_max << "] range for tolerated frequency ratio" << endl;
	}
	else
	{
		abs_freq_active = 0;
		abs_freq_inactive = 0;
		freq_ratio_min = 0;
		freq_ratio_max = 0;
	}
//	this->min_frequency = min_frequency;
//	this->min_frequency_per_class = min_frequency_per_class;

	time (&time_last_message);
}

MineDistances::~MineDistances() {
}

void MineDistances::free_memory(int index)
{
	if(DEBUG_OUT)
		cerr << "free memory for fragment: " << index << "\n" << flush;

	vector<int> * occ = data->get_smarts_occurences(index);
	for (unsigned int i = 0; i < occ->size(); i++)
	{
		if(DEBUG_OUT)
			cerr << "  free molecule: " << i << "\n" << flush;
		data->free_memory_from_id(occ->at(i),index);
		if ( mol_distances.find(occ->at(i)) != mol_distances.end() )
			mol_distances[occ->at(i)]->free_memory( index );
	}
}

bool MineDistances::calcDistance(int index1, int index2, vector<int> * occurences, bool check_freq)
{
	if (current_pair % 25 == 0)
	{
		time_t time_now;
		time (&time_now);
		if	( difftime (time_now,time_last_message) > 3 )
		{
			double percent = (current_pair/(double)num_max_pairs) * 100.0;
			fprintf(stderr,"%5.1f%s %9ld/%9ld\n",percent,"%",current_pair,num_max_pairs);
			time (&time_last_message);
		}
	}

	if(DEBUG_OUT)
		cerr << "    calculating distance: " << index1 << " <-> " << index2 << "\n" << flush;

	//OBSmartsPattern smartsPattern;

	vector<int> occurence_indices;
	vector<vector<double> > occurence_distances;

	// iterate over common smarts
	for (unsigned int i = 0; i < occurences->size(); i++)
	{
		//DEBUG_OUT = index1 == 2 && index2 == 16 && i == 51;

		/*
		OBMol * mol = data->get_mol( occurences->at(i) );

		smartsPattern.Init(*data->get_smarts(index1));
		smartsPattern.Match( *mol );
		vector<vector <int> > map1 = smartsPattern.GetUMapList();

		smartsPattern.Init(*data->get_smarts(index2));
		smartsPattern.Match( *mol );
		vector<vector <int> > map2 = smartsPattern.GetUMapList();
		 */

		vector<vector <int> > * map1 = data->get_matches_from_id(occurences->at(i), index1);//, true);
		vector<vector <int> > * map2 = data->get_matches_from_id(occurences->at(i), index2);

		if (map1->size() < 1 || map2->size() <1)
		{
			cerr << "\nERROR while calculating the distance: no match" << endl;
			cerr << "smiles: "<< *data->get_smiles_from_id(occurences->at(i)) << endl;
			cerr << "smarts 1: '" << *data->get_smarts(index1)<< "'" << endl;
			cerr << "matches: "<< map1->size() << endl;
			cerr << "smarts 2: '" << *data->get_smarts(index2)<< "'" << endl;
			cerr << "matches: "<< map2->size() << endl;

			OBMol * mol = new OBMol();
			OBConversion obconversion;
			obconversion.SetInFormat("smiles");
			obconversion.ReadString(mol,*data->get_smiles_from_id(occurences->at(i)));
			OBSmartsPattern smartsPattern;
			smartsPattern.Init( *data->get_smarts(index2) );
			smartsPattern.Match( *mol );
			vector<vector <int> > map = smartsPattern.GetMapList();
			cerr << "map-size of smarts 2: " << map.size() << endl;

			continue;
			//exit(1);
		}

		if(DEBUG_OUT)
		{
			cerr << "    mol " << (i+1) << "/" << occurences->size()<<"\n";
			cerr << "    frag 1 matches ";
			for (unsigned int j = 0; j < map1->size(); j++)
			{
				if (j>0)
					cerr << ", ";
				for (unsigned int k = 0; k < map1->at(j).size(); k++)
				{
					if (k>0)
						cerr << "-";
					cerr << map1->at(j)[k];
				}
			}
			cerr << endl;

			cerr << "    frag 2 matches ";
			for (unsigned int j = 0; j < map2->size(); j++)
			{
				if (j > 0)
					cerr << ", ";
				for (unsigned int k = 0; k < map2->at(j).size(); k++)
				{
					if (k>0)
						cerr << "-";
					cerr << map2->at(j)[k];
				}
			}
			cerr << endl;
		}

		if (DEBUG_OUT)
			cerr << "    mol: "<< *data->get_smiles_from_id(occurences->at(i)) <<"\n";

		if ( mol_distances.find(occurences->at(i)) == mol_distances.end() )
		{
			if (data->is_3d_enabled())
				mol_distances[occurences->at(i)] = new MolDistance3D( data->get_mol_from_id(occurences->at(i)) );
			else
				mol_distances[occurences->at(i)] = new MolDistance2D( data->get_mol_from_id(occurences->at(i)) );
		}

		vector<double> distances = mol_distances[occurences->at(i)]->get_distances( index1, map1, index2, map2 );

		if (distances.size()>0)
		{
			occurence_indices.push_back(occurences->at(i));
			occurence_distances.push_back(distances);
		}

		if (distances.size() >= 50 )
		{
			int unique[distances.size()];
			for (unsigned int w = 0; w < distances.size(); w++)
				unique[w] = -1;
			int uCount = 0;
			for (unsigned int w = 0; w < distances.size(); w++)
			{
				if (unique[w]==-1)
				{
					unique[w] = uCount;
					for (unsigned int y = w+1; y < distances.size(); y++)
						if (distances.at(w) == distances.at(y))
							unique[y] = uCount;

					uCount++;
				}
			}
			cerr << "WARNING: '" << *data->get_smarts(index1) << "' ("<< index1<<") and '" << *data->get_smarts(index2)
			<< "' ("<< index2<<") have " << distances.size() << " occurence-pairs (unique: "<< uCount << ") in '" << *data->get_smiles_from_id(occurences->at(i)) << "' ("<< occurences->at(i)<<")" << endl;

//			for (unsigned int w = 0; w < distances.size(); w++)
//				cerr << distances.at(w)<<"\n";
//			exit(1);
		}

		if(!CACHE_FRAGMENTS)
		{
			free_memory(index1);
			free_memory(index2);
		}
	}

	//cout << "X\n";

	if (check_freq && !checkMinFrequency(&occurence_indices, true))
		return false;

	cout << "[" << index1 << ";" << index2 << "]";
	if (occurence_indices.size() > 0)
		cout << ",{";
	for (unsigned int i = 0; i < occurence_indices.size(); i++)
	{
		if (i!=0)
			cout << ",";
		cout << "(" << occurence_indices[i] << ";";
		for (unsigned int j = 0; j < occurence_distances[i].size(); j++)
		{
			if (j!=0)
				cout << "#";
			if ((int)occurence_distances[i][j] == occurence_distances[i][j])
				cout << (int)occurence_distances[i][j];
			else
				cout << occurence_distances[i][j];
		}
		cout << ")";
	}
	if (occurence_indices.size() > 0)
		cout << "}";
	cout << endl;

	return true;
}


bool MineDistances::endsWithCChain(string * smarts)
{
	string chain ("C-C-C-C");

	bool endsWithCChain = false;
	if (smarts->size() > chain.size())
	{
		if (smarts->find(chain,smarts->size()-chain.size()) != string::npos)
			endsWithCChain = true;
		else
			endsWithCChain = (int)(smarts->find(chain))==0;
	}
	if (DEBUG_OUT && endsWithCChain)
		cerr << "    ends/start with "<< chain <<": "<< *smarts<<"\n";
	return endsWithCChain;
}


bool MineDistances::hasOnlyC(int smartIndex)
{
//	OBMol mol;
	OBSmartsPattern sp;

	// if (smarts_mols.find(smartIndex) == smarts_mols.end())
	// {

//		OBConversion obconversion;
//		obconversion.SetInFormat("smiles");
//		obconversion.ReadString(&mol, *data->get_smarts(smartIndex));
	sp.Init(*data->get_smarts(smartIndex));

	//	smarts_mols[smartIndex] = mol;
	//}
	//else
	//	mol = smarts_mols[smartIndex];


//	OBAtom atom;
//	FOR_ATOMS_OF_MOL(atom, mol)
//		if (!atom->IsCarbon())
//			return false;
    for(unsigned int i = 0; i < sp.NumAtoms(); i++)
    	//cerr << "atom " << i << " " << *data->get_smarts(smartIndex) << " is " << sp.GetAtomicNum(i) << "\n";
    	if (sp.GetAtomicNum(i) != 6)
    		return false;


//	if (DEBUG_OUT)
//		cerr << "    only C: "<<*data->get_smarts(smartIndex)<<"\n";
	return true;
}

bool MineDistances::checkMinFrequency(vector<int> * occurences, bool checkRatio)
{
	if (occurences->size() < (abs_freq_active + abs_freq_inactive))
	//if (occurences->size() < min_frequency)
	{
		if (DEBUG_OUT)
			cerr << "    not frequent: total " << occurences->size() << "/" << (abs_freq_active + abs_freq_inactive) << "\n";
		return false;
	}

	unsigned int activeCount = 0;
	for (unsigned int i = 0; i < occurences->size(); i++) {
		if ( data->is_active(occurences->at(i)))
			activeCount++;
	}
	unsigned int inactiveCount = occurences->size()-activeCount;

	if (activeCount < abs_freq_active || inactiveCount < abs_freq_inactive)
	//if (activeCount < min_frequency_per_class || inactiveCount < min_frequency_per_class)
	{
		if (DEBUG_OUT)
			cerr << "    not frequent: active " << activeCount << " / " << abs_freq_active << ", inactive " <<
				inactiveCount << " / " << abs_freq_inactive << "\n";
		return false;
	}

	if (checkRatio)
	{
		double ratio = activeCount / (double) inactiveCount;
		if (ratio < freq_ratio_min || ratio > freq_ratio_max)
		{
			if (DEBUG_OUT)
				cerr << "    frequency ratio violated: " << ratio << " not in [" << freq_ratio_min <<
					" - " << freq_ratio_max << "]" << "\n";
			return false;
		}
	}
	return true;
}

bool MineDistances::checkSubSmarts(string * smarts1, string * smarts2)
{
	if (smarts1->size() > smarts2->size())
		return smarts1->find(*smarts2) == string::npos;
	else
		return smarts2->find(*smarts1) == string::npos;
}

void MineDistances::mine( bool mineOneFragmentDistances )
{
	bool fragmentInvalid[data->num_smarts()];
	unsigned long validCount = 0;
	cerr << "Pre-check.. ";
	for (unsigned int i = 0; i < data->num_smarts(); i++)
	{
		if(hasOnlyC(i) || endsWithCChain(data->get_smarts(i)) || !checkMinFrequency(data->get_smarts_occurences(i),false))
			fragmentInvalid[i] = true;
		else
		{
			fragmentInvalid[i] = false;
			validCount++;
		}
	}
	cerr << "done ("<< validCount << "/" << data->num_smarts() << " fragments valid)" << endl;


	if (mineOneFragmentDistances)
		num_max_pairs = (validCount * (validCount + 1)) / (unsigned  long) 2;
	else
		num_max_pairs = (validCount * (validCount - 1)) / (unsigned  long) 2;
	current_pair = 0;

	unsigned long validPairCount = 0;

	for (unsigned int i = 0; i < data->num_smarts(); i++) {

		if (fragmentInvalid[i])
			continue;

		if(DEBUG_OUT)
  		  cerr << i << ": " << *data->get_smarts(i) << "\n";

		for (unsigned int j = i; j < data->num_smarts(); j++) {

			if (fragmentInvalid[j])
				continue;
			if (i==j && !mineOneFragmentDistances)
				continue;

			current_pair++;

			if(DEBUG_OUT)
 			  cerr << "  " << j << ": " << *data->get_smarts(j) << "\n";

			if (i!=j && !checkSubSmarts(data->get_smarts(i),data->get_smarts(j)))
			{
				if(DEBUG_OUT)
					cerr << "    "<< *data->get_smarts(i) << " and "<< *data->get_smarts(j) << " are sub smarts\n";
				continue;
			}
			if(DEBUG_OUT)
				cerr << "    joining fragments\n";

			vector<int> joined_occurences = join_occurences(i, j);

			if (!checkMinFrequency(&joined_occurences, true))
				continue;

			if (calcDistance(i,j,&joined_occurences, true))
				validPairCount++;
		}

		free_memory(i);
	}
	cerr << "Mining done ("<< validPairCount << "/" << num_max_pairs << " pairs valid)" << endl;
}

vector<int> MineDistances::join_occurences(int smarts_index_1, int smarts_index_2)
{
	vector<int> joined_occurences;
	vector<int> * occ_1 = data->get_smarts_occurences(smarts_index_1);
	vector<int> * occ_2 = data->get_smarts_occurences(smarts_index_2);

	for (unsigned int k = 0; k < occ_1->size(); k++)
	{
		bool match = false;
		for (unsigned int l = 0; l < occ_2->size(); l++)
		{
			if (occ_1->at(k) == occ_2->at(l))
			{
				match = true;
				break;
			}
		}
		if (match)
			joined_occurences.push_back(occ_1->at(k));
	}

	return joined_occurences;
}

void MineDistances::check()
{
	num_max_pairs = data->num_distance_pairs();

	for (current_pair = 0; current_pair < data->num_distance_pairs(); current_pair++)
	{
		if(DEBUG_OUT)
			cerr << "\nchecking pair: " << current_pair << " (of " << data->num_distance_pairs() << ")\n" << flush;
		int * idx = data->get_distance_pair(current_pair);
		if(DEBUG_OUT)
			cerr << "pair indices: " << idx[0] << " and " << idx[1] << " (of " << data->num_smarts() << ")\n" << flush;
		vector<int> joined_occurences = join_occurences(idx[0], idx[1]);
		calcDistance(idx[0],idx[1],&joined_occurences, false);
	}
}
