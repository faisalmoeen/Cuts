/*
#include "stdafx.h"
#include "gridline.h"


gridline::gridline(void)
{
}


gridline::~gridline(void)
{
}
*/
#include "stdafx.h"
#include "gridline.h"

inline int GridLine::leftBelowOf( int cid )
{
	int l = leftOf(cid);
	if(l>-1)
	{
		int lb = belowOf(l);
		return (lb>-1) ? lb : l;
	}	
	int b = belowOf(cid);
	return (b>-1) ? b : cid;
}
inline int GridLine::rightUpperOf( int cid )
{
	int r = rightOf(cid);
	if(r>-1)
	{
		int ru = upperOf(r);
		return (ru>-1) ? ru : r;
	}	
	int u = upperOf(cid);
	return (u>-1) ? u : cid;
}

/*---------------------------------------------------------------------
   DESCRIPTION: read the buffer and build a grid index
   AUTHOR: Hoyoung Jeung, 22/10/2007
   NOTE: typedef vector<MPoint*> GridNode;
		vector<GridNode> _grid; // 2*2 matrix
		typedef vector<MPoint> Trajectory;
		typedef map<int,Trajectory> Database; // <id, trajectory>
*---------------------------------------------------------------------*/
void GridLine::build( vector<PolyLine>* lns, bool star, double e, double tol)
{
	_lns = lns;
	_star = star;

	// compute the extent of the data space
	double maxarea =0;
	for(int i=0; i<_lns->size(); i++)
	{
		_ext.unionWith(_lns->at(i).ext);
		//maxarea = MAX(maxarea,_lns->at(i).ext.area());  // important for speeding up
	}

	//_csize = sqrt(maxarea);
	_csize = e + 2*tol;
	//_csize = e;

	// setting the grid size
	_ncol = ceil(_ext.width()/_csize);  
	_nrow = ceil(_ext.height()/_csize);

	if(_ncol*_nrow <4)
		return;

	if(_nodes.size()>0)
		_nodes.clear();

	for(int i=0; i<_nrow; i++)
		for(int j=0; j<_ncol; j++)
		{
			int c = i*_ncol + j;
			MBR m;
			cell2MBR(c,m);
			_mbrs[c] = m;
		}

	for(int lid=0; lid<_lns->size(); lid++)
	{
		set<int> cells;
		getCells(_lns->at(lid),0,cells);
		for(set<int>::iterator cell=cells.begin(); cell!=cells.end(); cell++)
			_nodes.insert(multimap<int,int>::value_type(*cell,lid));
	}
}

/*---------------------------------------------------------------------
    DESCRIPTION: compute all the cells intersecting given polyline.
				 actually results may contain some cells do not intersect
				 but no false negatives
    PARAMETERS 		
		1. p: given a point
    RETURN
		1. -1: if given point is out the extent 
    AUTHOR: Hoyoung Jeung, 6, Apr, 2007  
	NOTE: how to name cell ID? when row=2, col=5
	5 6 7 8 9 
	0 1 2 3 4
*---------------------------------------------------------------------*/
inline void GridLine::getCells( PolyLine& l, int bufsize, set<int>& rst )
{
	int c1 = getCell(&l.ext.min());
	int c2 = getCell(&l.ext.max());
	for(int i=1; i<bufsize; i++) // get more(bufsize) cells around the cells of given polyline
	{
		c1 = leftBelowOf(c1);
		c2 = rightUpperOf(c2);
	}
	int c1col = (c1<_ncol) ? c1 : c1%_ncol;
	int c2col = (c2<_ncol) ? c2 : c2%_ncol;
	int c1row = (c1-c1col)/_ncol;
	int c2row = (c2-c2col)/_ncol;
	int ncol = c2col - c1col + 1;
	int nrow = c2row - c1row + 1;

	// overlapping between the bounding box of l and cells
	int lloop = l.size()-1;
	for(int i=0; i<nrow; i++)
		for(int j=0; j<ncol; j++)
		{
			// removing deadspace from the bounding box overlapping areas
			int cellid = c1col+j+(c1row+i)*_ncol;
			//for(int z=0; z<lloop;z++)
			//	if(_mbrs[cellid].intersect(*l[z],*l[z+1]))
					rst.insert(cellid);
		}
}


inline void GridLine::getCells( MBR& m, int bufsize, set<int>& rst )
{
	int c1 = getCell(&m.min());
	int c2 = getCell(&m.max());
	for(int i=1; i<bufsize; i++) // get more(bufsize) cells around the cells of given polyline
	{
		c1 = leftBelowOf(c1);
		c2 = rightUpperOf(c2);
	}
	int c1col = (c1<_ncol) ? c1 : c1%_ncol;
	int c2col = (c2<_ncol) ? c2 : c2%_ncol;
	int c1row = (c1-c1col)/_ncol;
	int c2row = (c2-c2col)/_ncol;
	int ncol = c2col - c1col + 1;
	int nrow = c2row - c1row + 1;

	// overlapping between the bounding box of l and cells
	for(int i=0; i<nrow; i++)
		for(int j=0; j<ncol; j++)
			rst.insert(c1col+j+(c1row+i)*_ncol);
}

