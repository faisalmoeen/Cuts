#pragma once

#ifndef __GRIDLINE__
#define __GRIDLINE__

#include "grid.h"

class GridLine : public Grid
{
public:
	GridLine() { _ncol = _nrow = -1; };
	virtual ~GridLine() {};

	void getCells(PolyLine& l, int bufsize, set<int>& rst);
	void getCells(MBR& m, int bufsize, set<int>& rst);

	void	build(vector<PolyLine>* _lns, bool star, double e, double tol);
	void	rangeQuery(int qry, double r, set<int>& rst)
	{
		PolyLine* q = &(_lns->at(qry)); // query polyline
		set<int> cands; // candidate cells that overlap to the query range
						//int bufsize = 1 + int(r/_csize);
		getCells(*q, 1, cands);

		//double prune = 1 - (double)cands.size()/_lns->size();

		for (set<int>::iterator it = cands.begin(); it != cands.end(); it++) // for each candidate cell
		{
			// find all entries (points/pid) having the key (cell id)
			pair<GridEntry, GridEntry> entries = _nodes.equal_range(*it);

			for (GridEntry entry = entries.first; entry != entries.second; entry++)
			{
				int lid = (*entry).second;
				PolyLine* l = &(_lns->at(lid));
				double r2 = (r + q->tol + l->tol) * (r + q->tol + l->tol);
				double d2;
				q->ext.min_dist2(&l->ext, &d2);
				if (d2 <= r2)
				{
					if (_star)  // CuTS*
						q->distLL2_ST(l, &d2);
					else // CuTS
						q->distLL2(l, &d2);

					if (d2 <= r2)
						rst.insert(lid);
				}
			}
		}
	}

	int		rightUpperOf(int cid);
	int		leftBelowOf(int cid);

protected:

	vector<PolyLine>* _lns;
	bool			_star;
};


#endif 
