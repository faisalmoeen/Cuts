#pragma once

#ifndef __GRID__
#define __GRID__

#include "geom.h"
#include "buf.h"

class Grid
{
public:
	Grid() { _ncol = _nrow = -1; };
	virtual ~Grid() {};

	typedef multimap<int, int>::iterator GridEntry; //< cell id, pid (or lid)>

													// compute the center coordinates of cid cell and put the result into the given point
													//void	getCenter(int cid, Point &rst); 
	int		getCell(Point p); // compute the cell where the given point belongs to
	MBR*	getExtent() { return &_ext; };
	void	getNeighbors(int cid, set<int>& rst);
	void	cell2MBR(int cell, MBR& rst);

	void	build(vector<IPoint>* _pts, double e);
	void	rangeQuery(vector<IPoint>* _pts, int center, double radius, set<int>& rst);
	int		numCells() { return _nrow*_ncol; };
	bool	isBuilt() { return (_nodes.size()>0); }

protected:

	multimap<int, int> _nodes;  // index entries, each entry is <cell id, pid>
	map<int, MBR> _mbrs;  // <cell id, mbr>
	int		_ncol;	// # of columns
	int		_nrow;	// # of rows
	double	_csize;	// the width/height of each cell 
	MBR		_ext;	// entire extent of the map

	int	rightOf(int cid); // return the right cell of given cell
	int	upperOf(int cid); // return the upper cell of given cell
	int	belowOf(int cid); // return the below cell of given cell
	int	leftOf(int cid); // return the left cell of given cell
};


#endif 
