/*
#pragma once
class cvq
{
public:
	cvq(void);
	~cvq(void);
};

*/

#ifndef __CONVOYQUERY__
#define __CONVOYQUERY__

#include <time.h>
#include "buf.h"
#include "dbscan.h"
#include "cv.h"


class ConvoyQuery  
{
public:
	ConvoyQuery();
	virtual ~ConvoyQuery(){};

	int		_fneg; // # of false negatives
	int		_fpos; 
	int		_fneg2; // # of false negatives
	int		_fpos2; 
	int		_exac; // # of exact answers
	int		_corr; // # of correct answers
			
	virtual clock_t discover(int m, int k)= 0;	
	virtual clock_t discover(int m, int k, double e)= 0;	
	virtual void	summary( ofstream& fout )=0;
	void	write(ofstream& fo);
	void	write(char* fname); 
	void	showProgress(bool mode){_showprogress=mode;};	
	void	setRange(double e) {_e=e;}
	void	useIndex(int mode){_index=mode;}
	void	validate( vector<Convoy>& ans);
	void	output( Convoy& c );
	vector<Convoy>* getResults(){return &_rst;}


protected:

	int		_m;
	int		_k;
	double	_e;
	Buffer*	_buf;	/* inputs */		
	clock_t	_elapsed_t;	/* elapsed time */
	int		_index;
	clock_t	_index_t;
	bool	_showprogress;
	vector<Convoy>	_rst;	/* output */  
	
	double	findEps(int m, int k);
	void	initialize();	

};
#endif 
