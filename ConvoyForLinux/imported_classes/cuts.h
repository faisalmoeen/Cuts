/*#pragma once
class cuts
{
public:
	cuts(void);
	~cuts(void);
};
*/

#ifndef __CuTS__
#define __CuTS__

#include "cvq.h"
#include "gridline.h"
#include <list>
#include <math.h>
#include "idist.h"

#define NOINDEX		0	
#define IDISTANCE	1
#define GRIDLINE	2

class CuTS : public ConvoyQuery  // Convoy using Trajectory Simplification
{
public:
	CuTS(){defaultSetting();};
	CuTS(Buffer *b){_buf=b;defaultSetting();};
	virtual ~CuTS();

	clock_t	discover(int m, int k);
	clock_t	discover(int m, int k, double e){setRange(e);return discover(m,k);}
	clock_t	simplify(double tol);
	//clock_t simplify(){_tol = compute_tol();return simplify(_tol);}

	// setting	
	void	lambda(int glen){if(glen>0)_lamda=glen;}
	void	star(bool mode){_star=mode;};
	void	plus(bool mode){_plus=mode;};
	void	batch(int mode){_batch=mode;}  // 0: no batch, 1-semi batch, 2-batch
	void	actualTol(bool mode){_acttol=mode;}
	int		reduction(){return (int)_reduc;};
	double	getTolerance(){return _tol;}
	void	summary( ofstream& fout );
	void	summary2( ofstream& fout );
	void	print();
	void	matlab( int nth, char* outf);
	void	dumpPartition(char* fn);
	double	compute_tol(double e);	
	int		computeLambda(int m, int k, double reduc);

	class PolyLines : public vector<PolyLine>
	{
	public:
		MBR		_ext;
		void	push_back(PolyLine& l){_ext.unionWith(l.ext);vector<PolyLine>::push_back(l);}
	};


protected:    
 
	map<int,PolyLines> _tpart; // <id,polylines> time partition
	map<int,GridLine> _gline;  
	map<int,iDistance> _idist;
	int			_lamda;  // the number of time of each time partition 
	double		_tol;	// tolerance for simplification
	bool		_star;	// spatio-temporal mode	cuts-star
	bool		_plus;  // optimization for DP performance, cuts-plus
	bool		_acttol;
	double		_reduc; // reduction ratio.
	double		_effec;
	int			_batch;
	clock_t		_simple_t;
	clock_t		_refine_t;  // elapsed time of refinement step
	clock_t		_filter_t;  // elapsed time of  step
	vector<Convoy>	_cands;	/* candidates */  
	list<Point *> _splits;
	

	void	dp(vector<Point>& v, int j, int k, vector<double>& tols);
	void	compute_tols(vector<Point>& v, int j, int k, set<double>& order);
	
	void	refine(int m, int k, double e);
	void	verify(const Convoy& cand,vector<Convoy>& rst, int m, int k, double e);
	void	candidate(Convoy& c);
	void	analyze(vector<Convoy>& newrst);
	void	analyze2(vector<Convoy>& newrst);
	clock_t buildIndex();
	void	defaultSetting();
	

	
};
#endif 



