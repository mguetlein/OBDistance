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

bool DEBUG_OUT = false; // set via -v
bool CACHE_FRAGMENTS = true;

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

void compute3d(char *smiles)
{
	cerr << "Computing 3d for smiles: "<< smiles << endl;
    OBMol * mol = new OBMol();
    OBConversion obconversion;
    obconversion.SetInFormat("smiles");
    obconversion.ReadString(mol,smiles);
    cerr << "formula "<< mol->GetFormula() << "\n";
    cerr << "has3D: " << mol->Has3D() << "\n";
    OBBuilder builder;
    builder.Build(*mol);
    cerr << "has3D: " << mol->Has3D() << "\n";
    delete(mol);
}

void readSDF()
{
    OBMol * molSmi = new OBMol();
    OBConversion obconversionSmi;
    obconversionSmi.SetInFormat("smiles");
    obconversionSmi.ReadString(molSmi,"O=C(C(C(C=C3)=CC=C3O)=CO2)C1=C2C=C(O)C=C1O");
    cerr << "formula "<< molSmi->GetFormula() << "\n";

	  OBConversion obconversion;
	  obconversion.SetInFormat("sdf");
	  OBMol mol;

	  bool notatend = obconversion.ReadFile(&mol,"/home/martin/data/tmp/corinna_plz/nctrer_3D.sdf");
	  while (notatend)
	  {
		std::cout << "formula " << mol.GetFormula() << " Molecular Weight: " << mol.GetMolWt() << " 3D? " << mol.Has3D() << std::endl;

		mol.Clear();
		notatend = obconversion.Read(&mol);
	  }
}


int main(int argc, char **argv) {

//	readSDF();
//	return 0;

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


//	test();
//	exit(0);

	int status = 0;
	int c;

	char* smi_file = 0;
	char* fragment_file = 0;
	char* activity_file = 0;
	char* dist_file = 0;
	char* sdf_file = 0;
//	char* sdf_smiles_file = 0;

	bool check_smarts = false;
	bool mine_distance = false;
	bool aromaticity = false;
	bool check_distance = false;
	bool freq_set = false;
	bool enable_3d = false;
	bool one_frag_distances = false;


	char* compute_3d_from_smiles = 0;

//	int min_freq = 4;
//	int min_freq_per_class = 2;
	double rel_min_freq_per_class = 0.01;
	double freq_ratio_tolerance = 0.25;

	while ((c = getopt(argc, argv, "s:f:c:pdam:i:t:xzovr:e:q:l:")) != -1) {
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
		case 'l':
			sdf_file = optarg;
			break;
//		case 'y':
//			sdf_smiles_file = optarg;
//			break;
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
		case 'q':
			compute_3d_from_smiles = optarg;
			break;
		case ':':
			status = 1;
			break;
		case '?':
			status = 1;
			break;
		}
	}

	if (compute_3d_from_smiles)
	{
		compute3d(compute_3d_from_smiles);
		return 0;
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
		if (!activity_file || !dist_file || freq_set)
			status = 1;
	}

	if (!enable_3d)
	{
		if (sdf_file)// || sdf_smiles_file)
			status = 1;
	}
//	if (sdf_file && !sdf_smiles_file)
//		status = 1;

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
		     << "\t-l <sdf file> (reads 3d structure of smiles-molecules from sdf file)\n"
//		     << "\t-y <sdf smiles file> (mandatory to assign sdf-molecules to smiles)\n"
		     << "\t-o (Mine one fragment distances) (default: off)\n"
		     << "\t-v (Verbose, i.e. print debug messages) (default: off)\n"
		     << "\t-q <smiles> (Debug option: compute 3d from smiles and exit) (default: off)\n"
		     << "\n"
		     << "examples:\n"
		     << " - Check fragments with aromaticity:\n"
		     << "   " << argv[0] << " -s file.smi -f file.linfrag -p -a\n"
			 << " - Mine distances:\n"
			 << "   " << argv[0] << " -s file.smi -f file.linfrag -c file.class -d -m 6 -i 3\n"
			 << " - Check distances:\n"
			 << "   " << argv[0] << " -s file.smi -f file.linfrag -c file.class -t file.dist -x\n"
			 << "\n"
			 << "extended options via environment variables:\n"
			 << " - disable fragment caching with OBD_CACHE_FRAGMENTS=0, default is 1\n";
		return (status);
	}

	Data * d = new Data(smi_file, fragment_file, aromaticity, enable_3d, activity_file, dist_file, sdf_file); //, sdf_smiles_file);
	cerr << "Aromaticity: "<< aromaticity << endl;
	cerr << "3D-Enabled: "<< enable_3d << endl;

	char * cache_fragments = getenv("OBD_CACHE_FRAGMENTS");
	if (cache_fragments!=NULL && strcmp(cache_fragments,"0")==0)
		CACHE_FRAGMENTS = false;
	cerr << "Fragment caching enabled: "<< CACHE_FRAGMENTS << endl;

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
