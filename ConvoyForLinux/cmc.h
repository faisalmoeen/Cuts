/*#pragma once
class cmc
{
public:
	cmc();
	~cmc();
};
*/

#pragma once

#ifndef __CMC__
#define __CMC__

#include "mc.h"

class CoherentMovingCluster : public MovingCluster
{
public:
	CoherentMovingCluster(Buffer *b);
	virtual ~CoherentMovingCluster() {};

	clock_t	discover(int m, int k);
	clock_t	discover(int m, int k, double e) {
		setRange(e);
		return discover(m, k);
	}

};
#endif
