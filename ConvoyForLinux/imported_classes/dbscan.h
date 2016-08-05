/*
#pragma once
class dbscan
{
public:
	dbscan(void);
	~dbscan(void);
};

*/

#ifndef __DBSCAN__
#define __DBSCAN__

#include "cv.h"
#include "grid.h"

#define UNCLASSIFIED (-2)
#define NOISE (-1)
#define CLASSIFIED_START (0)


/*
 *	Implementation of DBSCCAN algorithm which was introduced in
 *  SIGKDD96 paper.
 */

class Dbscan  
{
public:	
	Dbscan(vector<IPoint>* ptset, int minpts, double eps);			
	virtual ~Dbscan(){};
	
	/* perform dbscan algorithm */
	void	perform();
		
	void	setIndex(Grid* grid) {_grid = grid;};  // perform with the grid index

	/* return the number of clustered discovered */
	int		numClusters() {return _clsts.size();};

	/* return i-th cluster */
	Convoy* getCluster(int i) {return &_clsts.at(i);};

	/* return all the clusters discovered */
	vector<Convoy>* getClusterSet() {return &_clsts;};


protected:
	vector<IPoint>* _pts;	// input point set
	int			_minpts;	    // dbscan parameter
	double		_eps;	    // dbscan parameter
	vector<int>	_cids;	    // cluster id container
	vector<Convoy> _clsts;     // result clusters discovered
	Grid*		_grid;		// grid index
	
	/* main dbscan function */
	void	dbscan();

	/* expanding a cluster */
	bool	expandCluster(int Point, int cid);

	/* return a cid corresponding to oid */
	int		clId(int row) {return _cids[row];}

	/* assign cluster id value*/
	int		nextId(int cid) {return (cid==NOISE)?CLASSIFIED_START:++cid;}	

	/* set seed point with the given id value */
	void	changeClId(int row, int cid){_cids[row]=cid;}

	/* set all seed points with the given id value */
	void	changeClIds(set<int>& seeds, int cid);

	/* perform a range query with a given center point and radius */
	void regionQuery(int center, double radius, set<int>& rst);

	/* construct the result cluster set from cid array(cids) */
	void	buildClustSet();

	/* find the iterator (location) having given value */
	vector<int>::iterator find(vector<int> *v, int value);

};

#endif 
