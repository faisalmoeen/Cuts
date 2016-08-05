/*
#include "stdafx.h"
#include "cuts.h"


cuts::cuts(void)
{
}


cuts::~cuts(void)
{
}
*/

#include "stdafx.h"
#include "cuts.h"
#include <algorithm>
#include <cassert>
#include "func.h"
#include "dbscanline.h"
#include <stdexcept>
#include <iterator>

void CuTS::defaultSetting()
{	
	ConvoyQuery::initialize();
	_star = _plus = false; 
	_simple_t = _refine_t = _index_t = 0;
	_batch=1; 
	_reduc = 90.0;
	_lamda=-1;
	_tol = -1.0;
	_index=NOINDEX;
	_acttol=true;
}


CuTS::~CuTS()
{
	for(list<Point*>::iterator i=_splits.begin(); i!=_splits.end(); i++)
		delete *i;
}

/*---------------------------------------------------------------------
	DESCRIPTION: create grid indices for every timestamp.
	we build grid index as many as number timestamps
	- double e : e constraint
	- vector<,Grid> &_grids : empty container to store results
	AUTHOR: Hoyoung Jeung, 22/10/2007
*---------------------------------------------------------------------*/
clock_t CuTS::buildIndex()
{
	if(!_gline.empty())	_gline.clear();
	if(!_idist.empty())	_idist.clear();

	clock_t start = clock();

	map<int,PolyLines>::iterator mi= _tpart.begin();
	int i=0;
	while(mi!=_tpart.end())  // for each polyline groups, it is sorted by time
	{
		if(mi->second.size()>=4) 
		{
			if(_index==IDISTANCE)
				_idist[i].build(&mi->second,_star,_e,_tol); 
			else if(_index==GRIDLINE)
				_gline[i].build(&mi->second,_star,_e,_tol); 
		}
		mi++;
		i++;
	}
	return clock() - start;
}


/*---------------------------------------------------------------------
	DESCRIPTION: CuTS (Convoy over Simplified Trajectories)
		perform a convoy query over trajectory-simplification
		method. it performs line generalization first before
		doing DBSCAN
	PARAMETERS:
	1. st: spatio-temporal consideration for both line simplification
	and distance measure among simplified trajectories
	AUTHOR: Hoyoung Jeung, 11-20 June 2007
*---------------------------------------------------------------------*/
clock_t CuTS::discover(int m, int k)
{	
	initialize();	
	_cands.clear();
	_m = m; _k = k;	
	if(_e<0) _e = findEps(m,k);

	if(_index)
		_index_t = buildIndex();// ----- important
	
	_filter_t = clock();
	vector<Convoy> V; /* set of current clusters */	
	int n = _tpart.size();
	int s=0;	
	map<int,PolyLines>::iterator mi= _tpart.begin();
	int ltp = mi->first-1; // last time partition id
	while(s<n)  // for each polyline groups, it is sorted by time
	{

//		cout << "For each polyline group: " << s<< endl;
		vector<Convoy> Vnext;
		vector<Convoy>::iterator c,v;

		//cout << " * Pasa For each polyline group: " << s<< endl;
		
		if(mi->second.size()<m || ltp+1 != mi->first)
		{
		//cout << " * Enter if " << endl;

			for(v=V.begin(); v!=V.end(); v++)
			{
				//cout << "for(v=V.begin(); v!=V.end(); v++) /  v->lifetime() [" << v->lifetime() << "] >= k [" << k << "]";
				if(v->lifetime()>=k)	
					candidate(*v);					
			}
				V=Vnext;
			ltp=mi->first;
			mi++;
			s++;
			continue;
		}
		//cout<< "Begin DBscanline"<<endl;			
		DbscanLine ds(&(mi->second),m,_e,_star);  /*	DBSCAN for polylines */	
		//cout<< "End DBscanline"<<endl;
		
		if(_index==IDISTANCE)	ds.setIndex(&_idist[s]);	
		else if(_index==GRIDLINE)	ds.setIndex(&_gline[s]);	
		ds.perform();  

		vector<Convoy> L = *ds.getClusterSet();
		for(v=V.begin(); v!=V.end(); v++)	
		{
			//cout << "for(v=V.begin(); v!=V.end(); v++)	 "  <<endl;

			v->_assigned = false;
			for(c=L.begin(); c!=L.end(); c++)
			{	
				//cout << "c=L.begin(); c!=L.end(); c++ ";
				/* we may need consecutive time check here */
				vector<int> is(v->size());
				vector<int>::iterator ii=set_intersection(v->begin(),v->end(),
					c->begin(),c->end(), is.begin());				
				is.erase(ii, is.end()); /* trim rst to contain only results */		
				int is_size = is.size();			
				if(is_size>=m) /* verification of convoy */	
				{	
					//cout << "entra if is_size [" << is_size << "] > m[" << m <<"]";
					v->_assigned = true;					
					v->_te = (s==n-1) ? _buf->lasttime() : _buf->timeof((s+1)*(_lamda-1));
					//v->_te = mi->second._ext.maxt;
					v->clear();  /* Vnext = Vnext U (v intersect c) */
					v->insert(is.begin(),is.end());
					Vnext.push_back(*v);
					c->_assigned = true;	
				}   
			} 
			if(!v->_assigned)
			{
				if(v->lifetime()>=k)			
 					candidate(*v);
			}
		} 
		if(s==n-1)
		{
			for(v=V.begin(); v!=V.end(); v++)
				if(v->lifetime()>=k)	
					candidate(*v);	
			break;
		}

//	cout << "discover 2"<<endl;
		/* the entry of the active list */
		for(c=L.begin(); c!=L.end(); c++)
			if(!c->_assigned)	
			{
				c->_ts = _buf->timeof(s*(_lamda-1));
				c->_te = (s==n-1) ? _buf->lasttime() : _buf->timeof((s+1)*(_lamda-1));
				Vnext.push_back(*c);
			}

		V=Vnext; 
		ltp=mi->first;
		mi++;
		s++;
		if(_showprogress) Function::progress(s,n);
	}

	/* 3. Refinement Phase  */ 
	refine(m,k,_e);

	_filter_t = clock()-_filter_t;
	_elapsed_t = _filter_t + _simple_t + _refine_t;	
	analyze(_cands); // for debugging
	analyze2(_cands); // for debugging
	print(); 
	return _elapsed_t;
}

/*---------------------------------------------------------------------
    DESCRIPTION: compute the number of false negatives and false positives
    PARAMETERS  
        1. ans : a convoy set which is supposed to be
        2. cand : a candidate set for testing
        3. fp : result for the number of false positive
		4. fn : result for the number of false negative
    RETURNS
        1. # of the same convoys
    AUTHOR: Hoyoung Jeung, 20 June 2007
*---------------------------------------------------------------------*/
void CuTS::refine(int m, int k, double e)
{		
	if(_cands.size()<1) return;
	
	clock_t start = clock();
	vector<Convoy> newrst;
	_fpos = _fneg = _exac = _corr = 0;

	for(int i=0; i<_cands.size(); i++)  // for old results
	{
		vector<Convoy> rtn;
		verify(_cands[i],rtn,m,k,e);
		for(int z=0;z<rtn.size(); z++)
			_rst.push_back(rtn[z]);
	}		
	_refine_t= clock()-start;	
}


void CuTS::analyze(vector<Convoy>& cands)
{		
	if(cands.size()<1) return;
	_fpos = _fneg = _exac = _corr = 0;

	for(int i=0; i<_rst.size(); i++)  // for old results
	{
		for(int j=0; j<cands.size(); j++)
		{
			if(_rst[i]==cands[j])
				_exac++;
			else if(_rst[i]<=cands[j])	
				_corr++;
		}
	}  
	_fpos2 = (1.0 - (_exac+_corr)/(double)cands.size())*100;
	_fneg2 = (1.0 - (_exac+_corr)/(double)_rst.size())*100;
	_fpos = (1.0 - (_exac)/(double)cands.size())*100;
	_fneg = (1.0 - (_exac)/(double)_rst.size())*100;
}
void CuTS::analyze2(vector<Convoy>& cands)
{		
	_effec=0;
	for(int j=0; j<cands.size(); j++)
		_effec+= cands[j].size() * cands[j].size() * cands[j].lifetime();
}


/*---------------------------------------------------------------------
	DESCRIPTION: 
	PARAMETERS: 
		- vector<Convoy>& cand : 
		- vector<int>& fp : 
		- vector<int>& fn : 
		- vector<int>& ex : 
	AUTHOR: Hoyoung Jeung, 30/10/2007 
	RETURN
		1 : false positive
		0 : identical
		-1 : false negative
*---------------------------------------------------------------------*/
void CuTS::verify(const Convoy& cand, vector<Convoy>& rst, int m, int k, double e)
{
	int csize = cand.size();	
	vector<Convoy> V; /* set of current clusters */	
	int s=cand._ts;
	int be = _buf->beforeof(cand._te); // before ending
	while(s<=cand._te)
	{		
		vector<Convoy> Vnext;
		vector<Convoy>::iterator c,v;

		/*	DBSCAN for only objects that the candidate contains */	
		vector<IPoint>* tmp= _buf->objAtTime(s);		
		int n = tmp->size();
		if(n<m || s==be)
		{
			for(v=V.begin(); v!=V.end(); v++)
				if(v->lifetime()>=k)		
					output(*v);
			if(s==be) return;
			V=Vnext;
			s=_buf->nextof(s);
			continue;
		}

		vector<IPoint> pts; pts.reserve(csize);
		for(int i=0; i<n; i++)
			if(cand.find(tmp->at(i).o)!=cand.end())
				pts.push_back(tmp->at(i));
	
		Dbscan ds(&pts,m,e);
		ds.perform();

		vector<Convoy> L = *ds.getClusterSet();
		for(v=V.begin(); v!=V.end(); v++)	
		{
			v->_assigned = false;	
			for(c=L.begin(); c!=L.end(); c++)
			{		
				vector<int> is(v->size());
				vector<int>::iterator intsecit=set_intersection(v->begin(),v->end(),c->begin(), c->end(),is.begin());				
				is.erase(intsecit, is.end()); /* trim rst to contain only results */	
				int is_size = is.size();						
				if(is_size>=m)/* verification of convoy */				
				{					
					v->_assigned = true;
					v->_te = s; 
					v->clear(); /* Vnext = Vnext U (v intersect c) */
					v->insert(is.begin(), is.end());
					Vnext.push_back(*v);
					c->_assigned = true;					
				}
			}
			if(!v->_assigned && v->lifetime()>=k)					
				output(*v);	
		}

		for(c=L.begin(); c!=L.end(); c++)	
			if(!c->_assigned)		
			{
				c->_ts= c->_te = s;
				Vnext.push_back(*c);	
			}

		V=Vnext;	
		s=_buf->nextof(s);
	}
}



/*---------------------------------------------------------------------
   DESCRIPTION: compute a proper tolerance 
	  - set<double> order :
   AUTHOR: Hoyoung Jeung, 27/11/2007
   NOTE: 
*---------------------------------------------------------------------*/
double CuTS::compute_tol(double e)
{
	int ntrial = (int) _buf->num_obj()/10; // 10 %

	Database::iterator it = _buf->getDB()->begin();   
	double tolsum=.0;
	int cnt=0;
	for(int i=0;i<ntrial; i++) // check the first ntrial numbers of trajectories
	{
		set<double> sorted_tols; // storage for tolerances sorted
		compute_tols(it->second, 0, it->second.size()-1, sorted_tols);

		int n = sorted_tols.size();
		if(n<1)
		{
			it++;
			continue;
		}
		set<double>::iterator tol = sorted_tols.begin();
		double maxdv=-1;
		double rst, lasttol=*tol;
		for (int z=0; z<n; z++)
		{
			if (*tol>e)
				break;

			double dv = *tol - lasttol;
			//double dv = *tol/lasttol;
			if (maxdv<dv)
			{
				maxdv = dv;
				rst = *tol;
			}
			
			lasttol = *tol;
			tol++;
		}

		tolsum += rst;
		cnt++;
	}
	return tolsum/cnt;
}


int CuTS::computeLambda(int m, int k, double reduc)
{
	//double nm = _buf->num_tstamp() * _buf->num_obj() - _buf->num_points();
	//if(nm==0)
	//	return 2;
	//return 2 * _buf->num_tstamp() * _buf->num_obj() / nm * 100/reduc;

	//double nm = _buf->num_tstamp() *  _buf->num_obj() * _buf->num_obj()/1620;
	//if(nm==0)
	//	return 2;

	int c = (int)((1-k/_buf->avglen())*10);
	if(c<=0)
		c =1;
	return 2 * c * 100/reduc;
}

/*---------------------------------------------------------------------
	DESCRIPTION: This is the Douglas-Peucker recursive simplification 
				 routine. It just marks vertices that are part of the 
				 simplified polyline for approximating the polyline 
				 subchain v[j] to v[k].
    PARAMETERS  
		1. tol: approximation tolerance. if this value is high then
		        given line will be simpler.
	RETURN
	    1. elapsed time
    AUTHOR:
		1. Copyright 2002, softSurfer (www.softsurfer.com)
		   This code may be freely used and modified for any purpose
		   providing that this copyright notice is included with it.
	       SoftSurfer makes no warranty for this code, and cannot be held
           liable for any real or imagined damage resulting from its use.
           Users of this code must verify correctness for their application.
		2. modified by Hoyoung Jeung, 27-29 June 2006
*---------------------------------------------------------------------*/
clock_t CuTS::simplify(double tol)
{	
	if(_lamda<0)
	{
		//_plen = (int) _buf->num_tstamp()*(100-reduc_goal)/100; /** important */	
		//_plen = (int) _buf->num_tstamp()/(100-reduc_goal); /** important */
		_lamda=50;
	}
	if(!_tpart.empty()) 
		_tpart.clear();

	_simple_t= clock();	
	_tol = tol;
	//_tol = compute_tol();

	
	Database::iterator it = _buf->getDB()->begin();
	Database::iterator itend = _buf->getDB()->end();
	double ratiosum=0.0;	
		//cout << "Simplify 0 " << endl;
	//int ec = (int) _buf->num_obj()/1.5; // expected capacity (polylines) for each time partition.
	int ec = (int) _buf->num_obj(); // expected capacity (polylines) for each time partition.

	while ( it != itend )
	{
	//	cout << "Simplify main it" << endl;
		Trajectory* trj = &(it->second); // original trajectory
		int n = trj->size();			
		int oid = it->first;

		vector<double> tols(n); // tolerance recorder
		dp(it->second, 0, n-1, tols);
														 
		int npt=1; // # of simplified points
		Point* lp = &(trj->at(0)); // last point
		double ltol=0.0; // maximum tolerance of each polyline
		int tp1, tp2; 
		PolyLine l(oid); l.push_back(lp); 	
		tp1 = _buf->tstampof(lp->t)/(_lamda-1); // last time partition
		bool virgin=true;

		for (int j=1; j<n; j++) // for each point on an original trajectory
		{	
			//cout << "Simplify each point in original traj [ for j: "<< j << endl;
			if((tols[j]<=_tol) && (j!=n-1))
				ltol = MAX(ltol,tols[j]);
			else
			{
		//				cout << "Simplify 3 " << endl;
				Point* cp = &(trj->at(j)); // current position				
				tp2 = _buf->tstampof(cp->t)/(_lamda-1); // current time partition
				l.push_back(cp);

				if(tp1<tp2 || j==n-1)
				{					
					int tend = (j==n-1 && _buf->tstampof(cp->t) !=_buf->num_tstamp()-1 ) ? tp2+1:tp2;
					for(; tp1<tend; tp1++)
					{
		//cout << "Simplify 4 " << endl;
						if(_tpart[tp1].capacity()<ec)
							_tpart[tp1].reserve(ec);

						int ps = tp1*(_lamda-1); // partition start timestamp
						int pe = (tp1+1)*(_lamda-1); // partition end timestamp						

						PolyLine l_(oid);
						int lb = 0; // lower bound
						int ub = l.size()-1; // upper bound
						for(int i=0; i<l.size();i++)
						{ 
							//cout << "Simplify 5 " << endl;
							int t_ = _buf->tstampof(l[i]->t);
							if(virgin)
								continue;
							else if(t_<ps)
								lb = (i==0) ? 0 : MAX(lb,i);
							else if(ps==t_)
								lb = i;

							if(j==n-1)
								ub= l.size()-1;
							else if(pe<t_)
								ub = (i==l.size()-1) ? l.size()-1 : MIN(ub,i);
							else if(pe==t_)
								ub = i;
						}

						for(int k=lb; k<=ub;k++)					
							l_.push_back(l[k]);

						//if(_star)
						//{
						//	if(_buf->tstampof(l_[0]->t)<ps)  // split : need to be free
						//	{
						//		_splits.push_back(PolyLine::split(_buf->timeof(ps),l_[0],l_[1]));
						//		l_[0] = _splits.back();
						//		l_.computeExtent();
						//	}
						//	int end = l_.size()-1;
						//	if(_buf->tstampof(l_[end]->t)>pe)
						//	{
						//		_splits.push_back(PolyLine::split(_buf->timeof(pe),l_[end-1],l_[end]));
						//		l_[end] = _splits.back();
						//		l_.computeExtent();
						//	}
						//}
						if(_acttol)
						{
							l_.tol = ltol; 
							ltol = .0; // individual tolerance
						}
						else
							l_.tol = _tol;  // global tolerance
						_tpart[tp1].push_back(l_);
						virgin=false;
					}

			//				cout << "Simplify 7 " << endl;
					PolyLine ll(oid); 
					//ll.reserve(_plen); // mem alloc for speeding up
					l=ll; // maybe faster than clear() ;
					if(_buf->tstampof(lp->t)<(tp1*(_lamda-1)))
						l.push_back(lp);
					l.push_back(cp);					
				}

				lp = cp;								
				npt++;				
			}

		} 
		ratiosum += (double) npt/(double) n;
		it++;
//		cout << "Simplify it++ ratiosum: " << ratiosum << endl;
	}	
	_reduc = 100 - (double)(ratiosum/(double) _buf->num_obj()) * 100; // %	
	_simple_t = clock()-_simple_t;
	//f		cout << "Simplify fin " << endl;
    return _simple_t; 
}

/*---------------------------------------------------------------------
   DESCRIPTION: debugging purpose
   PARAMETERS
	  - char* fn : file name
   AUTHOR: Hoyoung Jeung, 15/12/2007
*---------------------------------------------------------------------*/
void CuTS::dumpPartition( char* fn)
{
	ofstream fout(fn, ios::out);
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(0);

	Database::iterator it = _buf->getDB()->begin();
	Database::iterator itend = _buf->getDB()->end();
	double ratiosum=0.0;	

	// writing simplified trajectory
	int  cnt=0;
	set<int> oids;
	fout << "simplified trajectory" << endl;
	while ( it != itend ) // for each object
	{
		Trajectory* trj = &(it->second); // original trajectory
		int n = trj->size();			
		oids.insert(it->first);

		vector<double> tols(n); // tolerance recorder for bottomup dp
		dp(it->second, 0, n-1, tols);

		//fout <<"oid " << it->first << " : " << trj->at(0).t ;
		fout <<"oid " << it->first << " : " << _buf->tstampof(trj->at(0).t) ;
		for (int j=1; j<n; j++) // for each point on an original trajectory
		{	
			if(tols[j]>_tol)
				//fout << " " << trj->at(j).t;
				fout << " " << _buf->tstampof(trj->at(j).t);
		}
		//fout <<" " << trj->at(n-1).t << endl;
		fout <<" " << _buf->tstampof(trj->at(n-1).t) << endl;
		cnt++;
		it++;
	}
	
	// writing partitions of the simplified trajectory
	fout << endl << endl;
	map<int,PolyLines>::iterator pi=_tpart.begin();
	Trajectory tr;
	fout << "partition length (# timestamps) " << _lamda << endl;
	while(pi!=_tpart.end())
	{
		int et =  (pi->first+1)*(_lamda-1);
		fout << "part " << pi->first << " : " << pi->second.size() << " lines, "
			 " timestamp[" << pi->first*(_lamda-1) << "," << (pi->first+1)*(_lamda-1) << "]" << endl;
		for(int i=0; i<pi->second.size();i++) // # of polylines
		{
			PolyLine* l = &(pi->second[i]);
			if(oids.find(l->o)==oids.end())
				continue;

			fout << "\t<oid " << l->o << ">";
			for(int j=0; j<l->size(); j++) // # of points in each polyline
				//fout <<" " << l->at(j)->t ;
				fout <<" " <<  _buf->tstampof(l->at(j)->t) ;
			fout << endl;			
		}
		pi++;
		fout << endl;
	}
	fout.close();
}



/*---------------------------------------------------------------------
	DESCRIPTION: This is the Douglas-Peucker recursive simplification 
				 routine. It just marks vertices that are part of the 
				 simplified polyline for approximating the polyline 
				 subchain v[j] to v[k].
				 When tc=false, it is the original DP algorithm, otherwise
				 it is spatio-temporal version of the DP algorithm 
				 , called TD-ST(top down spatiotemporal), which was 
				 introduced in "spatiotemporal compression techniques 
				 for moving point objects"-EDBT05 paper.
    PARAMETERS  
		1. v: polyline array of vertex points 
		2-3. j,k: indices for the subchain v[j] to v[k]
		4. mk[]: array of markers matching vertex array v[]	
    AUTHOR:
		1. Copyright 2002, softSurfer (www.softsurfer.com)
		   This code may be freely used and modified for any purpose
		   providing that this copyright notice is included with it.
	       SoftSurfer makes no warranty for this code, and cannot be held
           liable for any real or imagined damage resulting from its use.
           Users of this code must verify correctness for their application.
		2. modified by Hoyoung Jeung, 27-29 June 2006
*---------------------------------------------------------------------*/
inline void CuTS::compute_tols(vector<Point>& v, int j, int k, set<double>& order)
{
	if (k <= j+1) return;// there is nothing to simplify
	int     maxi = j;          // index of vertex farthest from S	
	double   maxd2 = 0;         // distance squared of farthest vertex
	int		ceni = j + (int) (k-j)/2; // index of the center of S. it's for dp-plus
	int		mincend = (int) (k-j)/2;  // index of the minimum distance to ceni
	int		minceni = -1;
	Point u,w,Pb;	
	u = v[k] - v[j]; 
	double cu = DOT(u,u);     // segment length squared	
	double  b, cw, dv2;        // dv2 = distance v->at(i) to S squared

	for (int i=j+1; i<k; i++)
	{
		if(_star) // spatio-temporal => time ratio distance from the beginning point 
		{
			b = (double)(v[i].t - v[j].t)/(double)(v[k].t - v[j].t); 
			Pb = v[j] + u * b; 
			dv2 = PDIST2(v[i], Pb);
		}
		else // spatial => original DP
		{		
			w = v[i] - v[j]; 
			cw = DOT(w,u);

			if ( cw <= 0 ) // current point is left of the segment
				dv2 = PDIST2(v[i], v[j]);
			else if ( cu <= cw ) // current point is right of the segment
				dv2 = PDIST2(v[i], v[k]);
			else { // current point is over the segment
				b = cw / cu;
				Pb = v[j] + u * b; 
				dv2 = PDIST2(v[i], Pb);
			}        
		}

		if(_plus)
		{
			int idist = ABS(ceni-i);
			if(idist < mincend)
			{
				mincend = idist;
				minceni = i;
			}
		}

		// test with current max distance squared
		if (dv2 <= maxd2) 
			continue;        

		maxi = i;
		maxd2 = dv2;
	} 
	
	if(maxd2>0)
	{
		if(_plus && minceni>-1)
			maxi=minceni;
		order.insert(sqrt(maxd2));		
		compute_tols( v, j, maxi, order );  // polyline v->at(j) to v->at(maxi)
		compute_tols( v, maxi, k, order );  // polyline v->at(maxi) to v->at(k)
	}
}


inline void CuTS::dp(vector<Point>& v, int j, int k, vector<double>& tols)
{
	if (k <= j+1) return;// there is nothing to simplify
	int     maxi = j;          // index of vertex farthest from S	
	double   maxd2 = 0;         // distance squared of farthest vertex
	int		ceni = j + (int) (k-j)/2; // index of the center of S. it's for dp-plus
	int		mincend = (int) (k-j)/2;  // index of the minimum distance to ceni
	int		minceni = -1;
	Point u,w,Pb;	
	u = v[k] - v[j]; 
	double cu = DOT(u,u);     // segment length squared	
	double  b, cw, dv2;        // dv2 = distance v->at(i) to S squared
	double tol2 = _tol*_tol;

	for (int i=j+1; i<k; i++)
	{
		if(_star) // spatio-temporal => time ratio distance from the beginning point 
		{
			b = (double)(v[i].t - v[j].t)/(double)(v[k].t - v[j].t); 
			Pb = v[j] + u * b; 
			dv2 = PDIST2(v[i], Pb);
		}
		else // spatial => original DP
		{		
			w = v[i] - v[j]; 
			cw = DOT(w,u);

			if ( cw <= 0 ) // current point is left of the segment
				dv2 = PDIST2(v[i], v[j]);
			else if ( cu <= cw ) // current point is right of the segment
				dv2 = PDIST2(v[i], v[k]);
			else { // current point is over the segment
				b = cw / cu;
				Pb = v[j] + u * b; 
				dv2 = PDIST2(v[i], Pb);
			}        
		}

		if(_plus)
		{
			int idist = ABS(ceni-i);
			if(idist < mincend)
			{
				mincend = idist;
				minceni = i;
			}
		}

		// test with current max distance squared
		if (dv2 <= maxd2) 
			continue;        

		maxi = i;
		maxd2 = dv2;
	} 

	if(maxd2>tol2)
	{
		if(_plus && minceni>-1)
			maxi=minceni;		

		tols[maxi] = sqrt(maxd2);
		dp( v, j, maxi, tols );  // polyline v->at(j) to v->at(maxi)
		dp( v, maxi, k, tols );  // polyline v->at(maxi) to v->at(k)
	}
	else  
	{
		for (int i=j+1; i<k; i++)   // putting individual tolerances
			tols[i]=sqrt(maxd2);
	}
}

/*---------------------------------------------------------------------
	DESCRIPTION: print the query results
	AUTHOR: Hoyoung Jeung, 25/10/2007
*---------------------------------------------------------------------*/
void CuTS::print()
{
	cout.setf(ios_base::fixed,ios_base::floatfield);
	cout.precision(1);
	double filter = (_elapsed_t==0) ? 0 : (double)_filter_t/_elapsed_t*100;
	double simple = (_elapsed_t==0) ? 0 : (double)_simple_t/_elapsed_t*100;
	double refine = (_elapsed_t==0) ? 0 : (double)_refine_t/_elapsed_t*100;
	cout.precision(1);
	cout << "CuTS";
	if(_star) cout << "*"; 
	if(_plus) cout << "+"; 
	cout << " (m="<<_m<<" k="<<_k<<" e="<<_e<< " tol " << _tol << ") : " 
		<< _rst.size() << "/" << _cands.size() << " convoys (==" << _corr << "  -" << _fneg << "),  " 		
		<< _elapsed_t/1000 << " sec (simple "<<_simple_t/1000 << ", idx (" << _index << ") "
		<< _index_t/1000 << "), reduc " 
		<< _reduc << "%, part_len " << _lamda 
		//<< "  (filter " << filter <<"%, refine " << refine << "%, simpl " << simple 
		//<< _tpart.size()  << " group (len=" << _glen << ")" 
		<< endl;
}

/*---------------------------------------------------------------------
   DESCRIPTION: write the results to a filestream
	  - ofstream& fout : file stream to write
   AUTHOR: Hoyoung Jeung, 25/10/2007
*---------------------------------------------------------------------*/
void CuTS::summary2( ofstream& fout )
{
	fout.setf(ios_base::fixed,ios_base::floatfield);
	fout.precision(1);
	//double cand = (_ncand==0) ? 0 : (double)_rst.size()/_ncand*100;
	double elapsed = (_elapsed_t==0) ? 0 : _elapsed_t/1000;
	double filter = (_filter_t==0) ? 0 : _filter_t/1000;
	double simple = (_simple_t==0) ? 0 : _simple_t/1000;
	double refine = (_refine_t==0) ? 0 : _refine_t/1000;
	fout << "CuTS";
	if(_star) fout << "*"; 
	if(_plus) fout << "+"; 
	fout << " (m="<<_m<<" k="<<_k<<" e="<<_e<< " tol " << _tol << ") : "; 
	fout<< "\t";
	fout << _rst.size() << "/" << _cands.size() << " convoys, " << elapsed 
		 << " sec (filter " << filter <<", refine " << refine << ", simpl " << simple 
		 << ", " << _reduc << "% reduc, " << "index " << _index_t << " sec, " 
		 << _tpart.size()  << " parts (lambda=" << _lamda << ")"; 		 
	if(_fneg>0 || _fpos>0 || _exac >0 || _corr>0)
		fout << ", exact " << _exac << ", correct " << _corr << ", f.pos " << _fpos << "%,"<< " f.neg " << _fneg << "% ";
	if(_effec)
		fout << ", refinement cost " << _effec;
	fout <<endl;
}

void CuTS::summary( ofstream& fout )
{
	fout.setf(ios_base::fixed,ios_base::floatfield);
	fout.precision(1);
	//double cand = (_ncand==0) ? 0 : (double)_rst.size()/_ncand*100;
	double elapsed = (_elapsed_t==0) ? 0 : _elapsed_t/1000;
	double filter = (_filter_t==0) ? 0 : _filter_t/1000;
	double simple = (_simple_t==0) ? 0 : _simple_t/1000;
	double refine = (_refine_t==0) ? 0 : _refine_t/1000;
	fout << "CuTS";
	if(_star) fout << "*"; 
	if(_plus) fout << "+"; 
	fout << " (m="<<_m<<" k="<<_k<<" e="<<_e<< " tol " << _tol << ") : "; 
	fout<< "\t";
	fout << _rst.size() << "/" << _cands.size() << " convoys, " << elapsed 
		<< " sec (filter " << filter <<", refine " << refine << ", simpl " << simple 
		<< ", " << _reduc << "% reduc, " << "index " << _index_t << " sec, " 
		<< _tpart.size()  << " parts (lambda=" << _lamda << ")"; 		 
	if(_fneg>0 || _fpos>0 || _exac >0 || _corr>0)
		fout << ", exact " << _exac << ", correct " << _corr << ", f.pos " << _fpos << "%,"<< " f.neg " << _fneg << "% ";
	if(_effec)
		fout << ", refinement cost " << _effec;

	/* showing cadidate details */
	fout << "candidates " << endl;
	for(int i=0; i<_cands.size(); i++)
		//_cands[i].write(fout);
	{
		fout << "k=" << _cands[i]._te-_cands[i]._ts+1 << "[" 
			<< _cands[i]._ts << "," <<_cands[i]._te <<"](["
			<< _buf->tstampof(_cands[i]._ts) << "," <<_buf->tstampof(_cands[i]._te) <<"]), m="
			<< _cands[i].size() << "[";		
		for(set<int>::iterator it = _cands[i].begin(); it != _cands[i].end(); it++)
		{	
			if(it != _cands[i].begin())
				fout << ",";
			fout << *it ;
		}	
		fout << "]" << endl;
	}
	fout << endl ;
	fout << "results " << endl;
}

/*---------------------------------------------------------------------
	DESCRIPTION: 
	PARAMETERS: 
		- Convoy& c : 
	RETURNS:
        - void
	AUTHOR: Hoyoung Jeung, 12/11/2007 
	NOTE: 
*---------------------------------------------------------------------*/
void CuTS::candidate( Convoy& c )
{
	if(_cands.empty())
	{
		_cands.push_back(c);
		return;
	}
	
	if(_batch==0)        // no batch
		_cands.push_back(c);
	else if(_batch==1)   // semi-batch processing
	{ 
		for(int j=0; j<_cands.size(); j++)
		{
			if( (c._te < _cands[j]._ts) || (_cands[j]._te < c._ts))
				continue;   

			set<int> is;
			set_intersection(_cands[j].begin(),_cands[j].end(), c.begin(),c.end(),
				insert_iterator<set<int>>(is,is.begin()));	
			//if (is.size()==c.size())
			int threshold = (int)(c.size()*0.8);
			if (is.size()>=threshold)  // 90% overlapping
			{
				_cands[j]._ts = MIN(_cands[j]._ts,c._ts);
				_cands[j]._te = MAX(_cands[j]._te,c._te);						
				return;
			}
		}
		_cands.push_back(c);
	}	
	else if(_batch==2)  // batch
	{
		for(int j=0; j<_cands.size(); j++)
		{
			if( (c._te < _cands[j]._ts) || (_cands[j]._te < c._ts))
				continue;   

			_cands[j] += c;
			return;
		}
		_cands.push_back(c);
	}	
}

/*---------------------------------------------------------------------
	DESCRIPTION: create matlab command to see the effect of simplifications
	PARAMETERS: 
		- int nobj : # of objects to see
		- int ntstamp : # of timestamps to see
	AUTHOR: Hoyoung Jeung, 24/11/2007 
*---------------------------------------------------------------------*/
void CuTS::matlab( int nth, char* outf)
{
	vector<Trajectory> trjs;

	int cnt=0;
	int oid=-1;
	Database::iterator di = _buf->getDB()->begin();
	while(di!= _buf->getDB()->end())
	{
		if(cnt++==nth)
		{
			oid = di->first;
			trjs.push_back(di->second);
			break;
		}
		di++;
	}

	map<int,PolyLines>::iterator pi=_tpart.begin();
	Trajectory tr;

	while(pi!=_tpart.end())
	{
		for(int i=0; i<pi->second.size();i++)
		{
			PolyLine* l = &(pi->second[i]);
			if(l->o == oid)
			{
				for(int j=0; j<l->size(); j++)
				{
					Point p = *(l->at(j));
					if(tr.size()==0)
						tr.push_back(p);
					else if (p==tr[tr.size()-1])
						continue;
						
					tr.push_back(p);
					
				}
				break;
			}
		}
		pi++;
	}
	trjs.push_back(tr);
	Buffer::matlab(trjs,outf);
}
