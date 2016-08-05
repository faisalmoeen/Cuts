/*#include "stdafx.h"
#include "mc.h"


mc::mc(void)
{
}


mc::~mc(void)
{
}
*/

#include "stdafx.h"
#include "mc.h"
#include <algorithm>
#include "dbscan.h"
#include <cassert>
#include "func.h"


MovingCluster::MovingCluster(Buffer *b)
{
	_buf=b;
	initialize();
}

/*---------------------------------------------------------------------
    DESCRIPTION: perform a convoy query with the MC2 algorithm which
				 was introduced in "On discovering moving clusters
				 in spatio-temporal data"-SSTD05 paper. 
	PARAMETERS  
		1. theta: theta value (A intersect B)/(A unionWith B)  
    AUTHOR: Hoyoung Jeung, 22-23 June 2006, (modified 4 June 2007)
	NOTE: the original pseudo code of MC2 in the paper has bugs at 
	      while(k>0) clause. we modified it.
*---------------------------------------------------------------------*/
clock_t MovingCluster::mc2( double theta,int m, int k)
{	
	initialize();
	_timedb_t = _buf->_timedb_t;

	_elapsed_t = clock();	

	if(_e<0) _e = findEps(m,k);
	if(_index && _index_t==0)
		_index_t = buildIndex(_e);
	_m = m; _k = k;
	
	vector<Convoy> G; /* set of current clusters */	
		
	int n = _buf->num_tstamp();
	for(int s=0; s<n; s++)
	{
		vector<Convoy> Gnext;
		vector<Convoy>::iterator c,g;

		// less number of points in this timestamp
		vector<IPoint>* pts= _buf->objAtTstamp(s);
		if(pts->size()<m)
		{
			for(g=G.begin(); g!=G.end(); g++)
				//if(g->lifetime()>=k)			
				if(g->lifetime()>=2)
					output(*g);
			G=Gnext;
			continue;
		}

		/*	DBSCAN */
		Dbscan ds(pts,m,_e);
		if(_index) ds.setIndex(&_grids[s]);	
		ds.perform();

		/* adding a candidate set to result */
		vector<Convoy> L = *ds.getClusterSet();
		
		for(c=L.begin(); c!=L.end(); c++)		
			c->_assigned = false;	

		for(g=G.begin(); g!=G.end(); g++)	
		{
			g->_assigned = false;

			int trial=0; int cid=-1;		

			/* pick a random object and check if it belongs to c */
			while(g->size()-trial>=m && cid==-1)
			{
				/* pick a random object (first point) from a cluster */		
				Convoy::iterator tmp=g->begin();
				for(int z=0;z<trial;z++)
					tmp++;				
				int oid=*tmp;
	
				/* find a cluster in L contains the object */							
				int counter=0;
				for(c=L.begin(); c!=L.end(); c++)
				{					
					if(!(c->find(oid)==c->end()))
					{
						cid = counter;
						break;
					}
					counter++;
				}
				trial++;
			}			
			
			if(cid>-1) /* when a cluster found */				
			{												
				vector<int> intsec(g->size());
				vector<int> unio(g->size()+c->size());				
				vector<int> diff(g->size());		
				vector<int>::iterator intsecit, unioit, diffit;
				intsecit=set_intersection(g->begin(),g->end(),c->begin(), c->end(), intsec.begin());
				unioit=set_union(g->begin(),g->end(),c->begin(), c->end(), unio.begin());
				diffit=set_difference(g->begin(),g->end(),c->begin(), c->end(), diff.begin());
				
				/* trim rst to contain only results */
				intsec.erase(intsecit, intsec.end());
				unio.erase(unioit, unio.end());
				diff.erase(diffit, diff.end());
											
				/* verification of moving cluster */
				double ratio = (double)intsec.size()/(double)unio.size();
				if(ratio >= theta)	
				{					
					g->_assigned = true;
					g->_te = _buf->timeof(s);

					/* Gnext = Gnext U g o c, replace g to c  */
					g->clear();
					for(Convoy::iterator cc=c->begin(); cc!=c->end(); cc++)
						g->insert(*cc);	

					Gnext.push_back(*g);
					c->_assigned = true;	
				}

			} /* if */
				
			if(!g->_assigned)
				//if(g->lifetime()>=k)
				if(g->lifetime()>=2)
					output(*g);	
					
		} /* for (g) */
		
		for(c=L.begin(); c!=L.end(); c++)	
			if(!c->_assigned)		
			{
				c->_ts=c->_te=_buf->timeof(s);
				Gnext.push_back(*c);
			}
			
		
        if(s==n-1)
        {
            for(g=Gnext.begin(); g!=Gnext.end(); g++)
                if(g->_assigned && g->lifetime()>=2)		
				//if(g->_assigned && g->lifetime()>=k)
                      output(*g);
            break;
        }
		
		G=Gnext;
		if(_showprogress) Function::progress(s,n);			
	}
	
	_elapsed_t = clock()-_elapsed_t + _timedb_t;
	cout << "MC2 (theta=" << theta << ") - "; print();
	return _elapsed_t;
}


/*---------------------------------------------------------------------
	DESCRIPTION: create grid indices for every timestamp.
	we build grid index as many as number timestamps
	- double e : e constraint
	- vector<,Grid> &_grids : empty container to store results
	AUTHOR: Hoyoung Jeung, 22/10/2007
*---------------------------------------------------------------------*/
clock_t MovingCluster::buildIndex(double e)
{	
	if(e==0)
		return 0;

	int n = _buf->num_tstamp();
	if(!_grids.empty())
		_grids.clear();

	clock_t start = clock();

	for(int i=0; i<n; i++)	
		_grids[i].build(_buf->objAtTstamp(i),e);

	return clock() - start;
}



/*---------------------------------------------------------------------
	DESCRIPTION: 
	PARAMETERS: 
		- ofstream& fout : 
	AUTHOR: Hoyoung Jeung, 7/11/2007 
	NOTE: 
*---------------------------------------------------------------------*/
void MovingCluster::summary( ofstream& fout )
{
	fout << " (m="<<_m<<" k="<<_k<<" e="<<_e << ") : "; 
	fout << _rst.size() << " convoys, " << _elapsed_t/1000 <<  " sec (timedb " << _timedb_t/1000;
	if(_index)	
		fout << ", index (" << _index << ") " << _index_t/1000 << " sec";
	if(_fneg>0 || _fpos>0 || _exac >0 || _corr>0)
		fout << ", exact " << _exac << ", correct " << _corr << ", f.pos " << _fpos << "%,"
		<< " f.neg " << _fneg << "% ";
	fout << ")" << endl;
}


void MovingCluster::print()
{
	cout << "(m="<<_m<<" k="<<_k<<" e="<<_e<< ") : " << _rst.size() << " convoys, " 
		<< _elapsed_t/1000 << " sec (timedb " << _timedb_t/1000 << ", idx (" << _index << ") "
		<< _index_t/1000 << ")\n";
	
}