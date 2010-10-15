/*
 * MineDistances.h
 *
 *  Created on: 28.07.2009
 *      Author: martin
 */

#ifndef MINEDISTANCES_H_
//#include <time.h>
//#include "CheckSmarts.h"
#include "MolDistance.h"
#include "Data.h"
//#include "OBDistance.h"
#define MINEDISTANCES_H_

using namespace std;
using namespace OpenBabel;



class MineDistances {
public:
	MineDistances(Data * data,
			double rel_min_freq_per_class=0.01, double freq_ratio_tolerance=0.25);
			//unsigned int min_frequency=4, unsigned int min_frequency_per_class=2);

	virtual ~MineDistances();

	void mine( bool mineOneFragmentDistances );
	void check();

private:
	Data * data;

	unsigned int abs_freq_active;
	unsigned int abs_freq_inactive;
	double freq_ratio_min;
	double freq_ratio_max;

//	unsigned int min_frequency;
//	unsigned int min_frequency_per_class;

	time_t time_last_message;
	unsigned long num_max_pairs;
	unsigned long current_pair;

	map<int, MolDistance *> mol_distances;

	vector<int> join_occurences(int smarts_index_1, int smarts_index_2);
	bool hasOnlyC(int smartIndex);
	bool endsWithCChain(string * smarts);
	bool checkMinFrequency(vector<int> * occurences, bool checkRatio);
	bool checkSubSmarts(string * smarts1, string * smarts2);
	bool calcDistance(int index1, int index2, vector<int> * occurences, bool check_freq);
};

#endif /* MINEDISTANCES_H_ */
