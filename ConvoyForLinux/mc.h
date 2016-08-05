/*
#pragma once
class mc
{
public:
	mc();
	~mc();
};
*/
#pragma once


#ifndef __MOVINGCLUSTER__
#define __MOVINGCLUSTER__

#include "cvq.h"

class MovingCluster : public ConvoyQuery
{
public:
	MovingCluster() {};
	MovingCluster(Buffer *b);
	virtual ~MovingCluster() {};

	clock_t	discover(int m, int k) { return mc2(0.8, m, k); }
	clock_t	discover(int m, int k, double e) { setRange(e);return mc2(0.8, m, k); }
	clock_t	mc2(double theta, int m, int k);
	void	print();
	void	summary(ofstream& fout);


protected:

	map<int, Grid> _grids;
	clock_t	_timedb_t;

	clock_t	buildIndex(double e);



};
#endif 



