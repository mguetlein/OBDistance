/*
 * OBDistance.cpp
 *
 *  Created on: 28.07.2009
 *      Author: martin
 */



#include "CheckSmarts.h"
#include "MineDistances.h"
#include "Data.h"
#include "MolDistance.h"
#include "assert.h"

using namespace std;
using namespace OpenBabel;

bool DEBUG_OUT = false;

void test()
{
	vector<string> smiles;
	/*smiles.push_back("c1ccc(O)ccc1");
	smiles.push_back("C-C-O");
	smiles.push_back("C1=C-C=C(O)-C=C1");*/
	//smiles.push_back("c1ccccc1");
	//smiles.push_back("c1ccc(N)cc1");
	//smiles.push_back("C=C-C-N");
	//smiles.push_back("C1=C-C=C(N)-C=C1");
	//smiles.push_back("C1CCNCC1");

	smiles.push_back("C-[Pb]");
	//smiles.push_back("Cc1ccc(cc1)C3CCCc2cncn23");

	vector<string> smarts;
	/*smarts.push_back("C");
	smarts.push_back("N");
	smarts.push_back("C-N");
	smarts.push_back("c:N");
	smarts.push_back("C=C-N");
	smarts.push_back("c:c");
	smarts.push_back("c:c:N");
	smarts.push_back("C-C-N");
	smarts.push_back("C-C=C-N");
	smarts.push_back("C=C-C-N");*/
	/*smarts.push_back("C=C=C-N");
	smarts.push_back("C-C-C-N");*/

	smarts.push_back("[Pb]");
	smarts.push_back("C-N-C");
	smarts.push_back("C-N=C");
	smarts.push_back("C=N=C");
	smarts.push_back("CNC");
	smarts.push_back("C:N:C");

	for (unsigned int i = 0; i < smiles.size(); i++)
	{
		cout << "" << smiles[i] << endl;
		printf("%12s %3s %3s\n","aromaticity","ON","OFF");

		for (unsigned int k = 0; k < smarts.size(); k++)
		{
			bool match_on = false;
			bool match_off = false;

			{
				OBMol mol;
				OBConversion obconversion;
				obconversion.SetInFormat("smiles");
				obconversion.ReadString(&mol, smiles[i]);

				OBSmartsPattern smartsPattern;
				smartsPattern.Init( smarts.at(k) );
				match_on = smartsPattern.Match( mol );
			}
			{
				OBMol mol;
				OBConversion obconversion;
				obconversion.SetInFormat("smiles");
				obconversion.ReadString(&mol, smiles[i]);

				OBAtom a;
				FOR_ATOMS_OF_MOL(a,mol)
					a->UnsetAromatic();
				OBBond b;
				FOR_ATOMS_OF_MOL(b,mol)
					b->UnsetAromatic();
				mol.SetAromaticPerceived();

				OBSmartsPattern smartsPattern;
				smartsPattern.Init( smarts.at(k) );
				match_off = smartsPattern.Match( mol );
			}

			printf("%12s %3d %3d\n",smarts.at(k).c_str(),match_on,match_off);
		}
		cout << endl;
	}

}


int main(int argc, char **argv) {

//	OBMol * mol = new OBMol();
//	OBConversion obconversion;
//	obconversion.SetInFormat("smiles");
//	obconversion.ReadString(mol, "C1(N=C(NC(C)C)N=C(N=1)OC)NC(C)C");
//	cout << "formula "<< mol->GetFormula() << "\n";
//	OBSmartsPattern smartsPattern;
//	smartsPattern.Init( "N" );
//	smartsPattern.Match( *mol );
//	vector<vector <int> > map = smartsPattern.GetMapList();
//	cout << "num matches of N: "<< map.size() << "\n";

//	OBMol * mol = new OBMol();
//	OBConversion obconversion;
//	obconversion.SetInFormat("smiles");
//	//obconversion.ReadString(mol, "C1(N=C(NC(C)C)N=C(N=1)OC)NC(C)C");
//	obconversion.ReadString(mol, "OC1=CC3=C([C@@]2([H])CC[C@@]4(C)[C@](CC[C@@H]4O)([H])[C@@]([H])2[C@H](CCCCCCCCCS(CCCC(F)(F)C(F)(F)F)=O)C3)C=C1");
//	cout << "formula "<< mol->GetFormula() << "\n";
//	cout << "Has3D: " << mol->Has3D() << "\n";
//	OBBuilder builder;
//	builder.Build(*mol);
//	cout << "Has3D: " << mol->Has3D() << "\n";

	//test();
//	exit(0);

	int status = 0;
	int c;

	char* smi_file = 0;
	char* fragment_file = 0;
	char* activity_file = 0;
	char* dist_file = 0;

	bool check_smarts = false;
	bool mine_distance = false;
	bool aromaticity = false;
	bool check_distance = false;
	bool freq_set = false;
	bool enable_3d = false;
	bool one_frag_distances = false;

//	int min_freq = 4;
//	int min_freq_per_class = 2;
	double rel_min_freq_per_class = 0.01;
	double freq_ratio_tolerance = 0.25;

	while ((c = getopt(argc, argv, "s:f:c:pdam:i:t:xzovr:e:")) != -1) {
		switch (c) {
		case 's':
			smi_file = optarg;
			break;
		case 'f':
			fragment_file = optarg;
			break;
		case 'c':
			activity_file = optarg;
			break;
		case 't':
			dist_file = optarg;
			break;
		case 'p':
			check_smarts = true;
			break;
		case 'd':
			mine_distance = true;
			break;
		case 'x':
			check_distance = true;
			break;
		case 'a':
			aromaticity = true;
			break;
		case 'm':
			cerr << "min freq deprecated, use rel min frequency per class" << endl;
			exit(1);
//			min_freq = atoi(optarg);
//			freq_set = true;
			break;
		case 'i':
			cerr << "min freq per class deprecated, use rel min freq per class" << endl;
			exit(1);
//			min_freq_per_class = atoi(optarg);
//			freq_set = true;
//			break;
		case 'r':
			rel_min_freq_per_class = atof(optarg);
			freq_set = true;
			break;
		case 'e':
			freq_ratio_tolerance = atof(optarg);
			freq_set = true;
			break;
		case 'z':
			enable_3d = true;
			break;
		case 'o':
			one_frag_distances = true;
			break;
		case 'v':
			DEBUG_OUT = true;
			break;
		case ':':
			status = 1;
			break;
		case '?':
			status = 1;
			break;
		}
	}

	if (!smi_file | !fragment_file)
		status = 1;

	if (check_smarts)
	{
		if (mine_distance || check_distance || freq_set || one_frag_distances)
			status = 1;
	}
	else if (mine_distance)
	{
		if (check_distance || !activity_file)
		{
			status = 1;
		}
		if (!enable_3d && one_frag_distances)
		{
			status = 1;
		}
	}
	else if (check_distance)
	{
		if (!activity_file || !dist_file)
			status = 1;
	}

	// print usage and examples for incorrect input
	if (status) {
		cerr << "usage: " << argv[0] <<"\n"
		     << "\t-s <smiles file> (mandatory)\n"
		     << "\t-f <fragment file> (mandatory)\n"
		     << "\t-c <class file>\n"
		     << "\t-t <distance pair file>\n"
		     << "\t-p (check fragments)\n"
		     << "\t-d (mine distances) (requires class-file)\n"
		     << "\t-x (check distances) (requires class-file and distance-file)\n"
		     << "\t-a (aromaticity on) (default: off)\n"
		     //<< "\t-m <minimum frequency> (default: 4, only used for mine distances)\n"
		     //<< "\t-i <minimum frequency per class> (default: 2, only used for mine distances)\n"
		     << "\t-r <relative minimum frequency per class> (default: 0.01, only used for mine distances)\n"
		     << "\t-e <frequency ratio tolerance> (default: 0.25, only used for mine distances)\n"
		     << "\t-z (Enable 3d) (default: off)\n"
		     << "\t-o (Mine one fragment distances) (default: off)\n"
		     << "\t-v (Verbose, i.e. print debug messages) (default: off)\n"
		     << "\n"
		     << "examples:\n"
		     << " - Check fragments with aromaticity:\n"
		     << "   " << argv[0] << " -s file.smi -f file.linfrag -p -a\n"
			 << " - Mine distances:\n"
			 << "   " << argv[0] << " -s file.smi -f file.linfrag -c file.class -d -m 6 -i 3\n"
			 << " - Check distances:\n"
			 << "   " << argv[0] << " -s file.smi -f file.linfrag -c file.class -t file.dist -x\n";
		     ;
		return (status);
	}

	Data * d = new Data(smi_file, fragment_file, aromaticity, enable_3d, activity_file, dist_file);
	cerr << "Aromaticity: "<< aromaticity << endl;
	cerr << "3D-Enabled: "<< enable_3d << endl;

	if (check_smarts)
	{
		CheckSmarts * check = new CheckSmarts(d);
		cerr << "Checking smarts" << endl;
		check->checkSmarts();
		//check->validateSmarts();
	}
	else if (mine_distance)
	{
		MineDistances * dist = new MineDistances(d, rel_min_freq_per_class, freq_ratio_tolerance); //  min_freq, min_freq_per_class);
		//cerr << "Minimum frequency: "<< min_freq << endl;
		//cerr << "Minimum frequency per class: "<< min_freq_per_class << endl;
		cerr << "Mining distances" << endl;
		cerr << "Mine one fragment distances: "<< one_frag_distances << endl;
		dist->mine(one_frag_distances);
	}
	else if (check_distance)
	{
		MineDistances * dist = new MineDistances(d);
		cerr << "Checking distances" << endl;
		dist->check();
	}

	cerr << "Exiting program" << endl;
	exit(0);
//	new CheckFragments(smi_file, fragment_file, activity_file, check);
}
