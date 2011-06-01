/*
 * MolDistance.h
 *
 *  Created on: 22.07.2010
 *      Author: martin
 */

#ifndef MOLDISTANCE_H_
//#include "Data.h"
#include "OBDistance.h"
#define MOLDISTANCE_H_

using namespace std;
using namespace OpenBabel;


class Pos3D
{
public:
	Pos3D(double x, double y, double z);
	virtual ~Pos3D();

	int weight;
	double x;
	double y;
	double z;

	double dist( Pos3D * p );
	void merge( int w_this, Pos3D * p, int w );
};

class MolDistance
{
public:
	virtual vector<double> get_distances(int frag1_key, vector<vector <int> > * frag1, int frag2_key, vector<vector <int> > * frag2 ) = 0;
	virtual void free_memory(int frag_key) {};
};

class MolDistance2D : public MolDistance {
public:
	MolDistance2D(OBMol * mol);
	virtual ~MolDistance2D();

	vector<double> get_distances(int frag1_key, vector<vector <int> > * frag1, int frag2_key, vector<vector <int> > * frag2 );
	void free_memory(int frag_key);
	void print(OBMol * mol);

private:
	int **distance;
	int *separate_index;
};

class MolDistance3D : public MolDistance {
public:
	MolDistance3D(OBMol * mol);
	virtual ~MolDistance3D();

	vector<double> get_distances(int frag1_key, vector<vector <int> > * frag1, int frag2_key, vector<vector <int> > * frag2 );
	void free_memory(int frag_key);

private:
	int *separate_index;
	OBMol * mol;
	map<int, vector<Pos3D *> * > frag_positions;
	map<int, vector< int > * > frag_sep_index;
	map<int, double > frag_max_radius;

	void mine_3d_positions( int frag_key, vector<vector <int> > * frag );

};

#endif
