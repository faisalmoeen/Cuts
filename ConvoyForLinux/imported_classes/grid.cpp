/*
#include "stdafx.h"
#include "grid.h"


grid::grid(void)
{
}


grid::~grid(void)
{
}
*/
#include "stdafx.h"
#include "grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "func.h"
#include <iostream>
#include <cassert>



/*---------------------------------------------------------------------
    DESCRIPTION: compute the cell where the given point belongs to				
				 this is an hash function
    PARAMETERS 		
		1. p: given a point
    RETURN
		1. -1: if given point is out the extent 
    AUTHOR: Hoyoung Jeung, 6, Apr, 2007  
	NOTE: how to name cell ID? when nrow=2, ncol=5
	5 6 7 8 9 
	0 1 2 3 4
*---------------------------------------------------------------------*/
inline int Grid::getCell(Point *p)
{
	int col =  (p->x - _ext.minx)/_csize;
	int row =  (p->y - _ext.miny)/_csize;

	return  row*_ncol + col; 
}

inline void Grid::cell2MBR(int cell, MBR& rst)
{
	int col = (cell<_ncol) ? cell : cell%_ncol;
	int row = (cell-col)/_ncol;
	rst.minx = col*_csize + _ext.minx;
	rst.miny = row*_csize + _ext.miny;
	rst.maxx = rst.minx + _csize;
	rst.maxy = rst.miny + _csize;
}

/*---------------------------------------------------------------------
   DESCRIPTION: read the buffer and build a grid index
   AUTHOR: Hoyoung Jeung, 22/10/2007
   NOTE: typedef vector<MPoint*> GridNode;
		vector<GridNode> _grid; // 2*2 matrix
		typedef vector<MPoint> Trajectory;
		typedef map<int,Trajectory> Database; // <id, trajectory>
*---------------------------------------------------------------------*/
void Grid::build(vector<IPoint>* pts, double cellsize)
{
	if(pts->size()<4)
		return;
	
	// compute the extent of the data space
	for(int i=0; i<pts->size(); i++)
		_ext.unionWith(pts->at(i).p); 

	_csize = cellsize;

	// setting the grid size
	_ncol = ceil(_ext.width()/_csize);  
	_nrow = ceil(_ext.height()/_csize);

	if(_nodes.size()>0)
		_nodes.clear();

	for(int pid=0; pid<pts->size(); pid++)
	{
		int nodeid = getCell(pts->at(pid).p);
		_nodes.insert(multimap<int,int>::value_type(nodeid,pid));
		if(_mbrs.find(nodeid)==_mbrs.end())
		{
			MBR m;
			cell2MBR(nodeid,m);
			_mbrs.insert(map<int,MBR>::value_type(nodeid,m));
		}
	}
}

/*---------------------------------------------------------------------
   DESCRIPTION: perform range query based on the grid
	  - int qry : id (the order of pts) of the center point
	  - double radius : range (should be e constraint)
	  - vector<int>* rst : result set
   AUTHOR: Hoyoung Jeung, 22/10/2007
   NOTE: how to name cell ID? when row=2, col=5
   5 6 7 8 9 
   0 1 2 3 4
   typedef vector<MPoint*>	GridNode;
*---------------------------------------------------------------------*/
void Grid::rangeQuery( vector<IPoint>* pts, int qry,double radius,set<int>& rst)
{
	Point* q = pts->at(qry).p; // query point 
	int		qryCell = getCell(q);
	set<int> cands; // candidate cells that overlap to the query range
	cands.insert(qryCell);
	getNeighbors(qryCell,cands);

	/* refinement */
	double radius2 = radius * radius;
	for(set<int>::iterator i=cands.begin(); i!=cands.end(); i++) // for each candidate cell
	{
		// find all entries (points/pid) having the key (nodeid/cell id)
		pair<GridEntry, GridEntry> entries = _nodes.equal_range(*i);

		//for (GridEntry entry=entries.first; entry != entries.second; entry++)
		for (GridEntry entry=entries.first; entry != entries.second; ++entry)
		{
			int pid =  (*entry).second;
			if(PDIST2(*q,*(pts->at(pid).p)) <= radius2)
				rst.insert(pid);
		}
	}
}

/*---------------------------------------------------------------------
   DESCRIPTION: search neighbor cells
   PARAMETERS
	  - int cid :
	  - set<int>& rst :
   AUTHOR: Hoyoung Jeung, 16/1/2008
*---------------------------------------------------------------------*/
inline void Grid::getNeighbors( int cid, set<int>& rst )
{
	int r=rightOf(cid); if(r>-1) rst.insert(r);
	int l=leftOf(cid);  if(l>-1) rst.insert(l);
	int u=upperOf(cid); if(u>-1) rst.insert(u);
	int b=belowOf(cid); if(b>-1) rst.insert(b);
	if(r>-1 && u>-1) rst.insert(upperOf(r)); // right-upper
	if(r>-1 && b>-1) rst.insert(belowOf(r)); // right-below
	if(l>-1 && u>-1) rst.insert(upperOf(l)); // left-upper
	if(l>-1 && b>-1) rst.insert(belowOf(l)); // left-below
}

/*---------------------------------------------------------------------
   DESCRIPTION: compute the right/left/upper/below cell of given cell
	  - int cid : cell id (node id) of this grid
   RETURN
	  -1 if there is no right cell (given cell is the right most)
   AUTHOR: Hoyoung Jeung, 23/10/2007
*---------------------------------------------------------------------*/
inline int Grid::rightOf( int cid )
{
	int rst = cid + 1;
	return (rst%_ncol==0) ? -1 : rst;
}
inline int Grid::leftOf( int cid )
{
	int rst = cid - 1;
	return (cid%_ncol==0) ? -1 : rst;
}
inline int Grid::upperOf( int cid )
{
	int rst = cid + _ncol;
	return (rst>=_nrow*_ncol) ? -1 : rst;
}
inline int Grid::belowOf( int cid )
{
	int rst = cid - _ncol;
	return (rst<0) ? -1 : rst;
}
