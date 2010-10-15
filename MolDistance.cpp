/*
 * MolDistance.cpp
 *
 *  Created on: 22.07.2010
 *      Author: martin
 */

#include "MolDistance.h"

using namespace std;
using namespace OpenBabel;

extern bool DEBUG_OUT;
bool DIST_DEBUG_OUT = DEBUG_OUT && false;


// used to store the positions of fragments
Pos3D::Pos3D(double x, double y, double z)
{
	this->weight = 1;
	this->x = x;
	this->y = y;
	this->z = z;
}

Pos3D::~Pos3D() {
}

double Pos3D::dist( Pos3D * p )
{
	return sqrt( pow((x - p->x), 2) +  pow((y - p->y), 2) + pow((z - p->z), 2) );
}

// w_this is not this->weight, set this->weight automatically
void Pos3D::merge( int w_this, Pos3D * p, int w )
{
	x = (x * w_this + p->x * w) / (double)( w_this + w );
	y = (y * w_this + p->y * w) / (double)( w_this + w );
	z = (z * w_this + p->z * w) / (double)( w_this + w );
}

// to enanble cerr << point << endl;
ostream& operator << (ostream& os,  Pos3D * p)
{
	char buffer[256];
	sprintf(buffer,"%5.2f/%5.2f/%5.2f",p->x,p->y,p->z);
	os << buffer;
	return os;
}

// for computing 3d distance between two fragments in a molecule
// fragments can occur multiple times in a molecule
// several merge mechanisms to reduce number of distances (else could be > 200)
MolDistance3D::MolDistance3D(OBMol * mol)
{
	this->mol = mol;

	// build 3d coordinates
	if (DIST_DEBUG_OUT)
		cerr << "    build 3D ...";
	OBBuilder builder;
	builder.Build(*mol);
	if (DIST_DEBUG_OUT)
		cerr << "    done" << endl;

	// set seprate_index:
	// each atom gets an index, same index <-> same structure in the OBMol object
	separate_index = new int[mol->NumAtoms()+1];
	int sep_count = 1;
	vector<OBMol> separatedMolecules = mol->Separate();
	for (unsigned int i = 0; i < separatedMolecules.size(); i++)
	{
		OBMol m = separatedMolecules.at(i);
		OBAtom a;
		FOR_ATOMS_OF_MOL(a,m)
		{
			separate_index[sep_count]=i;
			sep_count++;
		}
	}
}

MolDistance3D::~MolDistance3D() {
	delete(separate_index);
}

// sets frag_positions array for fragment, i.e. a 3d-pos for each occurence for frag in the molecule
// merge positions if distance closer than (maximum) radius of the fragment
void MolDistance3D::mine_3d_positions( int frag_key, vector<vector <int> > * frag )
{
	if ( frag_key==-1 || frag_positions.find(frag_key) == frag_positions.end() )
	{
		if (DIST_DEBUG_OUT)
			cerr << "      positions:\n";

		vector< Pos3D * > * pos = new vector< Pos3D * >();
		vector< int > * sep_index = new vector< int >();
		double max_radius = 0;


		for (unsigned int k = 0; k < frag->size(); k++)
		{
			Pos3D * p = new Pos3D(0,0,0);
			for (unsigned int i = 0; i < frag->at(k).size(); i++)
			{
				// compute position by merging atom position sequentially
				OBAtom * a = mol->GetAtom(frag->at(k)[i]);
				Pos3D * pp = new Pos3D( a->GetX(), a->GetY(), a->GetZ() );
				p->merge( i, pp, 1 );
				delete(pp);
			}
			pos->push_back( p );
			sep_index->push_back( separate_index[frag->at(k)[0]] );
			if (DIST_DEBUG_OUT)
				cerr << "        " << p << " [" <<separate_index[frag->at(k)[0]] << "]\n";

			// compute max_radius, i.e. highest distance from each atom to center position
			double max_dist = 0;
			for (unsigned int i = 0; i < frag->at(k).size(); i++)
			{
				OBAtom * a = mol->GetAtom(frag->at(k)[i]);
				Pos3D * pp = new Pos3D( a->GetX(), a->GetY(), a->GetZ() );
				double d = p->dist( pp );
				if (max_dist < d)
					max_dist = d;
				delete(pp);
			}
			if (max_dist > max_radius)
				max_radius = max_dist;
		}

		// merge positions only if there is more than one position
		// and the fragement size is > 1
		if (frag->size()>1 && frag->at(0).size()>1)
		{
			bool merge = true;
			while ( merge )
			{
				// loop: if merged once, try again, recompute all distances (PENDING)

				if (DIST_DEBUG_OUT)
					cerr << "        measure merge distances (radius: "<<max_radius<<"):\n";
				merge = false;
				for (unsigned int k = 0; k < pos->size()-1; k++)
				{
					for (unsigned int l = k+1; l < pos->size(); l++)
					{
						if (sep_index->at(k) != sep_index->at(l))
						{
							// positions are separated, never merge!
							if (DIST_DEBUG_OUT)
								cerr << "          "<<k<<"-"<<l<<" : SEPERATED\n";
							continue;
						}
						double d = pos->at(k)->dist( pos->at(l) );
						if (DIST_DEBUG_OUT)
							cerr << "          "<<k<<"-"<<l<<" : "<<d<<"\n";
						if (d < max_radius)
						{
							if (DIST_DEBUG_OUT)
								cerr << "          merging: "<< k << "-" << l << "\n";
							pos->at(k)->merge( pos->at(k)->weight, pos->at(l), 1 );
							pos->at(k)->weight++;
							pos->erase(pos->begin()+l);
							sep_index->erase(sep_index->begin()+l);
							merge = true;
							break;
						}
					}
					if (merge)
						break;
				}
				if (DIST_DEBUG_OUT)
				{
					if (merge)
					{
						cerr << "          after merging\n";
						for (unsigned int k = 0; k < pos->size(); k++)
							cerr << "            " << pos->at(k) << endl;
					}
				}
			}
		}
		frag_positions[frag_key] = pos;
		frag_sep_index[frag_key] = sep_index;
		frag_max_radius[frag_key] = max_radius;
	}
}


// compute distance between two fragments
// distances are ommitted if they are too small
// distances are merged if they are very similar
vector< double > *MolDistance3D::get_distances( int frag1_key, vector<vector <int> > * frag1, int frag2_key, vector<vector <int> > * frag2 )
{
	if (DIST_DEBUG_OUT)
		cerr << "    frag 1\n";
	mine_3d_positions(frag1_key, frag1);
	if (DIST_DEBUG_OUT)
		cerr << "    frag 2\n";
	mine_3d_positions(frag2_key, frag2);

	vector< double > * distances = new vector< double >();

	// use min-max-radius (i.e. take smaller radius of both (max-)radius) for merging, and as min distance
	double min_max_radius = min(frag_max_radius[frag1_key],frag_max_radius[frag2_key]);

	if (DIST_DEBUG_OUT)
		cerr << "    distances (radius: "<<min_max_radius<<")\n";

	// if same fragment, compute distance only if number of occurences > 2
	if (frag1_key == frag2_key && frag_positions[frag1_key]->size() < 2)
	{
		if (DIST_DEBUG_OUT)
			cerr << "      none\n";
	}
	else
	{
		for (unsigned int k = 0; k < frag_positions[frag1_key]->size(); k++)
		{
			for (unsigned int l = 0; l < frag_positions[frag2_key]->size(); l++)
			{
				// if same fragment, compute distance only once for each pair
				if (frag1_key == frag2_key && l <= k)
					continue;

				// check if same structure, i.e. not separated
				if (frag_sep_index[frag1_key]->at(k)==frag_sep_index[frag2_key]->at(l))
				{
					double d = frag_positions[frag1_key]->at(k)->dist(frag_positions[frag2_key]->at(l));
					// if distance is smaller than the smaller radius of both fragments, ommit it (overlap)
					if (d >= min_max_radius)
					{
						if (DIST_DEBUG_OUT)
							cerr << "      "<<k<<"-"<<l<<": "<<d<<"\n";
						distances->push_back(d);
					}
					else if (DIST_DEBUG_OUT)
						cerr << "      "<<k<<"-"<<l<<": "<<d<<" DISTANCE TOO SMALL\n";
				}
				else if (DIST_DEBUG_OUT)
					cerr << "      "<<k<<"-"<<l<<": SEPARATED\n";
			}
		}
	}

	// merge distances that are very similar
	if (distances->size()>1)
	{
		if (DIST_DEBUG_OUT)
			cerr << "    sort distances\n";
		std::sort(distances->begin(), distances->end());

		// if one fragment is a single-atom fragment the min_max_radius is 0
		// in this case merge distance with 0.25A differnce
		if (min_max_radius == 0)
			min_max_radius = 0.25;
		if (DIST_DEBUG_OUT)
			cerr << "    merge distances (radius: "<< min_max_radius <<")\n";

		int merge_weight = 1;
		for (unsigned int j = 0; j < distances->size()-1; j++)
		{
			if ( abs(distances->at(j) - distances->at(j+1)) < min_max_radius )
			{
				double m = (merge_weight * distances->at(j) + distances->at(j+1)) / (double)(merge_weight+1);
				if (DIST_DEBUG_OUT)
					cerr << "      merge: "<< j << "-" << (j+1) << " " <<
						 distances->at(j) << "(*"<< merge_weight << ") + " << distances->at(j+1) << " -> " << m << "\n";
				distances->at(j) = m;
				merge_weight++;
				distances->erase(distances->begin()+(j+1));
				j--;
			}
			else
				merge_weight = 1;
		}
		if (DIST_DEBUG_OUT)
		{
			for (unsigned int j = 0; j < distances->size(); j++)
				cerr << "      "<<distances->at(j)<<"\n";
		}
	}
	if (DIST_DEBUG_OUT)
		cerr << "\n";

	if (!DIST_DEBUG_OUT && DEBUG_OUT)
	{
		for (unsigned int j = 0; j < distances->size(); j++)
			cerr << "      "<<distances->at(j)<<"\n";
		if (DIST_DEBUG_OUT)
			cerr << "\n";
	}

	return distances;
}



MolDistance2D::MolDistance2D(OBMol * mol)
{
	distance = new int* [mol->NumAtoms()+1];
	for (unsigned int i = 1; i <= mol->NumAtoms(); i++)
		distance[i] = new int[mol->NumAtoms()+1];

	for (unsigned int i = 1; i <= mol->NumAtoms(); i++)
		for (unsigned int j = 1; j <= mol->NumAtoms(); j++)
			distance[i][j] = -1;

	//start index with 1
	separate_index = new int[mol->NumAtoms()+1];
	int sep_count = 1;
	vector<OBMol> separatedMolecules = mol->Separate();
	for (unsigned int i = 0; i < separatedMolecules.size(); i++)
	{
		OBMol m = separatedMolecules.at(i);
		OBAtom a;
		FOR_ATOMS_OF_MOL(a,m)
		{
			separate_index[sep_count]=i;
			sep_count++;
		}
	}

	OBAtom a;
	FOR_ATOMS_OF_MOL(a,mol)
	{
		for(OBMolAtomBFSIter b(mol, a->GetIdx()); b; b++ )
		{
			if (a->GetIdx() != b->GetIdx() && separate_index[a->GetIdx()]==separate_index[b->GetIdx()] )
			{
				distance[a->GetIdx()][b->GetIdx()] = b.CurrentDepth()-1;
			}
		}
	}
}

MolDistance2D::~MolDistance2D() {
}

void MolDistance2D::print(OBMol * mol)
{
		cerr << "  ";
		for (unsigned int i = 1; i <= mol->NumAtoms(); i++)
				cerr << i<< "(" << mol->GetAtom(i)->GetType() << ")" << " ";
		cerr << endl;

		for (unsigned int i = 1; i <= mol->NumAtoms(); i++)
		{
			cerr << i <<"(" << mol->GetAtom(i)->GetType() << ")"<< " ";
			for (unsigned int j = 1; j <= mol->NumAtoms(); j++)
			{
				if (j>1)
					cerr << ",";
				printf("%2d", distance[i][j] ); //get_distance(i, j));
			}
			cerr << endl;
		}
}

vector< double > *MolDistance2D::get_distances( int frag1_key, vector<vector <int> > * frag1, int frag2_key, vector<vector <int> > * frag2 )
{
	vector< double > * distances = new vector< double >();
	//if (DIST_DEBUG_OUT)
	//map->print(data->get_mol(occurences->at(i)));

	if (DIST_DEBUG_OUT)
		cerr << "    ";

	for (unsigned int j = 0; j < frag1->size(); j++)
	{
		for (unsigned int k = 0; k < frag2->size(); k++)
		{
			int atom_pair_count = 0;
			int occ_dist = INT_MAX;
			bool no_path = false;
			bool overlap = false;

			for (unsigned int l = 0; l < frag1->at(j).size(); l++)
			{
				for (unsigned int m = 0; m < frag2->at(k).size(); m++)
				{
					if (frag1->at(j)[l] == frag2->at(k)[m])
					{
						overlap = true;
						break;
					}
					else
					{
						int dist = distance[frag1->at(j)[l]][frag2->at(k)[m]]; //map->get_distance(frag1->at(j)[l], frag2->at(k)[m]);
						if (dist == -1)
						{
							no_path = true;
							break;
						}
						else
						{
							if (dist < occ_dist)
								occ_dist = dist;
							//occ_dist = (occ_dist*atom_pair_count + dist)/(double)(atom_pair_count+1);
							atom_pair_count++;
						}
					}
				}
				if(overlap || no_path)
					break;
			}

			if (DIST_DEBUG_OUT)
			{
				if(overlap)
					cerr<< "overlap, ";
				else if(no_path)
					cerr<<"no-path, ";
				else
					cerr << occ_dist << ", ";
			}

			if (!overlap && !no_path)
			{
				if (occ_dist == INT_MAX)
				{
					cerr << "WTF" << endl;
					exit(1);
				}

				bool match = false;
				vector<double>::iterator it;
				for (it = distances->begin(); it < distances->end(); it++)
				{
					if (*it == occ_dist)
					{
						match = true;
						break;
					}
					else if (*it > occ_dist)
						break;
				}
				if (!match)
					distances->insert(it, occ_dist);


				/*
				bool match = false;
				int insert = 0;
				for (unsigned int w = 0; w < distances.size(); w++)
				{
					if (distances[w] == occ_dist)
					{
						match = true;
						break;
					}
					if (distances[w] > occ_dist)
					{
						insert = w;
						break;
					}
				}
				if (!match)
					distances.insert(distances.begin().insert,occ_dist);
					*/
			}
		}
	}
	if (DIST_DEBUG_OUT)
		cerr << endl;
	return distances;
}
