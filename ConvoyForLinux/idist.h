#pragma once
#ifndef __iDISTANCE__
#define __iDISTANCE__

#include "geom.h"
#include "buf.h"

class iDistance
{
public:
	iDistance() { _tol = 0; };
	virtual ~iDistance() {};
	bool	isBuilt() { return (_dists_sorted.size()>0); }
	void	build(vector<PolyLine>* lns, bool star, double e, double tol)
	{
		_tol = tol;
		_lns = lns;
		_star = star;

		int nline = _lns->size();
		MBR ext;
		int mint = _lns->at(0).startTime();
		int maxt = _lns->at(0).endTime();
		for (int i = 0; i<nline; i++)
		{
			mint = MIN(mint, _lns->at(i).startTime());
			maxt = MAX(maxt, _lns->at(i).endTime());
			ext.unionWith(_lns->at(i).ext);
		}

		_ref = ext.center();


		_dists.reserve(nline);
		for (int lid = 0; lid<nline; lid++)
		{
			double d1;
			//if(star) 
			//	_ref.distLL_ST(&_lns->at(lid),&d1);
			//else
			//_ref.distLL(&_lns->at(lid),&d1);
			Point pp = _lns->at(lid).ext.center();
			d1 = PDIST(_ref, _lns->at(lid).ext.center());
			_dists_sorted.insert(multimap<double, int>::value_type(d1, lid));
			_dists.push_back(d1);
			_diago.push_back(_lns->at(lid).ext.diagonal() / 2);
		}
	}
	void	rangeQuery(int qry, double r, set<int>& rst)
	{
		PolyLine* q = &(_lns->at(qry)); // query polyline
		double R = r + q->tol + _tol;

		//int ncand=0;
		double ub = _dists[qry] + _diago[qry];
		double lb = _dists[qry] - _diago[qry];

		for (multimap<double, int>::iterator dist = _dists_sorted.begin(); dist != _dists_sorted.end(); dist++)
		{
			int lid = (*dist).second;

			PolyLine* l = &(_lns->at(lid));

			double R = r + q->tol + l->tol;
			double R2 = R*R;
			double d2;
			q->ext.min_dist2(&l->ext, &d2);

			ub += R + _diago[lid];
			lb -= (R + _diago[lid]);

			if ((lb <= _dists[lid]) && (_dists[lid] <= ub))
			{
				//ncand++;
				if (d2 <= R2)
				{
					if (_star)
						q->distLL2_ST(l, &d2);
					else
						q->distLL2(l, &d2);

					if (d2 <= R2)
						rst.insert(lid);
				}
			}
		}
		//double prune = (1 - ncand/(double)_dists_sorted.size())*100;
	}


protected:

	Point	_ref; // reference polyline
	double		_tol;
	vector<PolyLine>* _lns;
	bool		_star;

	multimap<double, int> _dists_sorted;  // index entries, each entry is <distance2.lid>
	vector<double> _diago;  // index entries, each entry is <distance2.lid>
	vector<double> _dists;  // index entries, each entry is <distance2.lid>
};


#endif 
