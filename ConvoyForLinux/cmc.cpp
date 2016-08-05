/*
cmc::cmc()
{
}


cmc::~cmc()
{
}
*/
//#include "stdafx.h"
#include "cmc.h"
#include <algorithm>
#include "func.h"

CoherentMovingCluster::CoherentMovingCluster(Buffer *b)
{
	_buf = b;
	initialize();
}


/*---------------------------------------------------------------------
DESCRIPTION: Coherent Moving Cluster (CMC)
perform a convoy query with the MC2 algorithm which
was introduced in "On discovering moving clusters
in spatio-temporal data"-SSTD05 paper. However, we
modify the algorithm to be able to find not moving
clusters but cv.
In convoy discovery, the theta value must be 1, thus
while clause in the original algorithm is not needed
here.
PARAMETERS
- int m : minimum number of objects
- int k : minimum lifetime
- double e : minimum closeness
AUTHOR: Hoyoung Jeung, 4 June 2007
TODO: do we need the continuousness time check? can it be pruned by
g->assigned ?
*---------------------------------------------------------------------*/
clock_t CoherentMovingCluster::discover(int m, int k)
{
	initialize();
	_timedb_t = _buf->_timedb_t;

	if (_e<0) _e = findEps(m, k);
	if (_index)
		_index_t = buildIndex(_e);
	_m = m; _k = k;

	_elapsed_t = clock();
	int n = _buf->num_tstamp();

	vector<Convoy> V; /* set of current clusters */
	int s = 0; // sweep line of timestamps
	while (s<n)
	{
		cout << "current time s: " << s << " - n: " << n << endl;
		vector<Convoy> Vnext;
		vector<Convoy>::iterator c, v;

		// less number of points in this timestamp
		vector<IPoint>* pts = _buf->objAtTstamp(s);
		// WE - If less points than m, just skip. no sence of doing it! JUST GO OUT with continue
		if (pts->size()<m)
		{
			for (v = V.begin(); v != V.end(); v++)
				if (v->lifetime() >= k)
					output(*v);
			s++;
			V = Vnext;
			continue;
		}

		/*	DBSCAN */
		Dbscan ds(pts, m, _e);
		if (_index) ds.setIndex(&_grids[s]);
		ds.perform();

		/* adding a candidate set to result */
		vector<Convoy> L = *ds.getClusterSet();

		for (v = V.begin(); v != V.end(); v++)
		{
			v->_assigned = false;
			for (c = L.begin(); c != L.end(); c++)
			{
				vector<int> is(v->size());
				vector<int>::iterator intsecit = set_intersection(v->begin(), v->end(), c->begin(), c->end(), is.begin());
				is.erase(intsecit, is.end()); /* trim rst to contain only results */
				int is_size = is.size();
				if (is_size >= m)/* verification of convoy */
				{
					v->_assigned = true;
					v->_te = _buf->timeof(s);
					v->clear(); /* Vnext = Vnext U (v intersect c) */
					v->insert(is.begin(), is.end());
					Vnext.push_back(*v);
					c->_assigned = true;
				}
			}
			if (!v->_assigned)
			{
				if (v->lifetime() >= k)
					output(*v);
			}
		} /* for (v) */

		if (s == n - 1)
		{
			for (v = V.begin(); v != V.end(); v++)
				if (v->lifetime() >= k)
					output(*v);
			break;
		}

		for (c = L.begin(); c != L.end(); c++)
		{
			if (!c->_assigned)
			{
				c->_ts = c->_te = _buf->timeof(s);
				Vnext.push_back(*c);
			}
		}

		V = Vnext;
		s++;
		if (_showprogress) Function::progress(s, n);
	}

	_elapsed_t = clock() - _elapsed_t + _index_t + _timedb_t;

	cout << "CMC "; print();
	return _elapsed_t;
}



