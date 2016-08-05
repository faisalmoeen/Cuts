#pragma once

#ifndef __DBSCAN_LINE__
#define __DBSCAN_LINE__

#include "cv.h"
#include "func.h"
#include "gridline.h"
#include "idist.h"
#include <set>

#ifndef UNCLASSIFIED
#define UNCLASSIFIED (-2)
#define NOISE (-1)
#define CLASSIFIED_START (0)
#endif 


class DbscanLine
{
public:
	DbscanLine(vector<PolyLine>* lns, int minobj, double eps, bool st);
	virtual ~DbscanLine() {};

	/* perform dbscan algorithm */
	void	perform();

	void	setIndex(GridLine* gl) { _gl = gl; }  // perform with the grid index
	void	setIndex(iDistance* id) { _idx = id; }  // perform with the grid index

													/* return the number of clustered discovered */
	int		numClusters() { return _clsts.size(); };

	/* return i-th cluster */
	Convoy* getCluster(int i) { return &_clsts.at(i); };

	/* return all the clusters discovered */
	vector<Convoy>* getClusterSet() { return &_clsts; };


	// lid is a order of lines in _lns. i.e., _lns[lid] indicates a line segment
	// cid is a cluster id.
protected:
	vector<PolyLine>* _lns;	// input point set
	int			_minlns;		// this is not MinLines, but min. # of objects
	double		_eps;		    // dbscan parameter
	vector<int>	_cids;		    // cluster id container
	vector<Convoy> _clsts;      // result clusters discovered
	bool		_star;          // measuring distance with considering spatial only or
	GridLine*	_gl;		// gridline index
	iDistance*	_idx;		// iDistance index

							/* main dbscan function */
	void	dbscan();

	/* expanding a cluster */
	bool	expandCluster(int Point, int cid);

	/* return a cid corresponding to oid */
	int		clId(int lid) { return _cids[lid]; }

	/* assign cluster id value*/
	int		nextId(int cid) { return (cid == NOISE) ? CLASSIFIED_START : ++cid; }

	/* set seed point with the given id value */
	void	changeClId(int lid, int cid) { _cids[lid] = cid; }

	/* set all seed points with the given id value */
	void	changeClIds(set<int>& seeds, int cid);

	/* perform a range query with a given center point and radius */
	void	regionQuery(int query, double radius, set<int>& rst);

	/* construct the result cluster set from cid array(cids) */
	void	buildClustSet();

};

#endif 
#ifndef __DBSCAN_LINE__
#define __DBSCAN_LINE__

#include "../cv/cv.h"
#include "../common/func.h"
#include "../idx/gridline.h"
#include "../idx/idist.h"
#include <set>

#ifndef UNCLASSIFIED
#define UNCLASSIFIED (-2)
#define NOISE (-1)
#define CLASSIFIED_START (0)
#endif 


class DbscanLine
{
public:
	DbscanLine(vector<PolyLine>* lns, int minobj, double eps, bool st);
	virtual ~DbscanLine() {};

	/* perform dbscan algorithm */
	void	perform();

	void	setIndex(GridLine* gl) { _gl = gl; }  // perform with the grid index
	void	setIndex(iDistance* id) { _idx = id; }  // perform with the grid index

													/* return the number of clustered discovered */
	int		numClusters() { return _clsts.size(); };

	/* return i-th cluster */
	Convoy* getCluster(int i) { return &_clsts.at(i); };

	/* return all the clusters discovered */
	vector<Convoy>* getClusterSet() { return &_clsts; };


	// lid is a order of lines in _lns. i.e., _lns[lid] indicates a line segment
	// cid is a cluster id.
protected:
	vector<PolyLine>* _lns;	// input point set
	int			_minlns;		// this is not MinLines, but min. # of objects
	double		_eps;		    // dbscan parameter
	vector<int>	_cids;		    // cluster id container
	vector<Convoy> _clsts;      // result clusters discovered
	bool		_star;          // measuring distance with considering spatial only or
	GridLine*	_gl;		// gridline index
	iDistance*	_idx;		// iDistance index

							/* main dbscan function */
	void	dbscan();

	/* expanding a cluster */
	bool	expandCluster(int Point, int cid);

	/* return a cid corresponding to oid */
	int		clId(int lid) { return _cids[lid]; }

	/* assign cluster id value*/
	int		nextId(int cid) { return (cid == NOISE) ? CLASSIFIED_START : ++cid; }

	/* set seed point with the given id value */
	void	changeClId(int lid, int cid) { _cids[lid] = cid; }

	/* set all seed points with the given id value */
	void	changeClIds(set<int>& seeds, int cid);

	/* perform a range query with a given center point and radius */
	void	regionQuery(int query, double radius, set<int>& rst);

	/* construct the result cluster set from cid array(cids) */
	void	buildClustSet();

};

#endif 
