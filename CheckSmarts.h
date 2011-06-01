/*
 * CheckFragments.h
 *
 *  Created on: 27.07.2009
 *      Author: martin
 */

#ifndef CHECKFRAGMENTS_H_
#include "Data.h"
#define CHECKFRAGMENTS_H_

using namespace std;
using namespace OpenBabel;

class CheckSmarts {

public:
	CheckSmarts(Data * data);
	virtual ~CheckSmarts();

	void checkSmarts();
	void validateSmarts();

private:
	Data * data;

};

#endif /* CHECKFRAGMENTS_H_ */
