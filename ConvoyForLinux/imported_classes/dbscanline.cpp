/*
#include "stdafx.h"
#include "dbscanline.h"

dbscanline::dbscanline(void)
{
}

dbscanline::~dbscanline(void)
{
}
*/

#include "stdafx.h"
#include "dbscanline.h"
#include "func.h"
#include <cassert>

DbscanLine::DbscanLine(vector<PolyLine> *lines, int minlines, double eps, bool st)
{ 
	_lns=lines;
	_minlns=minlines; 
	_eps=eps;
	_star=st;
	_gl=NULL;	
	_idx=NULL;	
}


/*---------------------------------------------------------------------
	DESCRIPTION: perform dbscan algorithm
    AUTHOR: Hoyoung Jeung, 7 June 2006
*---------------------------------------------------------------------*/
void DbscanLine::perform()
{
	/* initializing clusterID values */
	_cids.resize(_lns->size(),UNCLASSIFIED);   
	//cout << "enter dbscan"<<endl;	
	dbscan();
	//cout << "pass dbscan"<<endl;
	buildClustSet();
}

/*---------------------------------------------------------------------
	DESCRIPTION: perform the dbscan algorithm
    AUTHOR: Hoyoung Jeung, 25 oct 2005 (modified 6 June 2006)
    NOTES: psuedo code from the KDD96 paper below
	DBSCAN (SetOfPoints, Eps, MinPts)
		// SetOfPoints is UNCLASSIFIED
		ClusterId := nextId(NOISE);
		FOR i FROM 1 TO SetOfPoints.size DO
			Point := SetOfPoints.get(i);
			IF Point.ClId = UNCLASSIFIED THEN
				IF ExpandCluster(SetOfPoints, Point,ClusterId, Eps, MinPts) THEN
					ClusterId := nextId(ClusterId)
				END IF
			END IF
		END FOR
	END; // DBSCAN
*---------------------------------------------------------------------*/
void DbscanLine::dbscan()
{		
	int ClusterId = nextId(NOISE);	// the first cluster id begins from 0

	int n=_lns->size();

	for (int i=0; i<n; i++)
	{
		//	cout << " dbscan 1 i["<<i<<"] / n ["<<n<< "]"<<endl;	

		if (_cids[i]==UNCLASSIFIED) 
		{
			if (expandCluster(i, ClusterId))
			{
				ClusterId = nextId(ClusterId);				
			}
		}
	}
}

/*---------------------------------------------------------------------
	DESCRIPTION: 
    PARAMETERS  
		1. Point: row in db(the order in the SetOfPoints list)
		2. ClId: cluster id
    RETURN 
        1. true: expandable
		2. false: unexpndable
    AUTHOR: Hoyoung Jeung, 25,oct, 2005
    NOTES: psuedo code from the KDD96 paper below
	ExpandCluster(SetOfPoints, Point, ClId, Eps,MinPts) : Boolean;
		seeds:=SetOfPoints.regionQuery(Point,Eps);
		IF seeds.size<MinPts THEN // no core point
			SetOfPoint.changeClId(Point,NOISE);
			RETURN False;
		ELSE // all points in seeds are density-reachable from Point
			SetOfPoints.changeClIds(seeds,ClId);
			seeds.delete(Point);
			WHILE seeds <> Empty DO
				currentP := seeds.first();
				result := SetOfPoints.regionQuery(currentP, Eps);
				IF result.size >= MinPts THEN
					FOR i FROM 1 TO result.size DO
						resultP := result.get(i);
						IF resultP.ClId	IN {UNCLASSIFIED, NOISE} THEN
							IF resultP.ClId = UNCLASSIFIED THEN
								seeds.append(resultP);
							END IF;
							SetOfPoints.changeClId(resultP,ClId);
						END IF; // UNCLASSIFIED or NOISE
					END FOR;
				END IF; // result.size >= MinPts
				seeds.delete(currentP);
			END WHILE; // seeds <> Empty
			RETURN True;
		END IF
	END; // ExpandCluster
*---------------------------------------------------------------------*/
inline bool DbscanLine::expandCluster(int Point, int ClId)
{
	set<int> seeds;
	regionQuery(Point, _eps, seeds);	

	//cout << " expandCluster 1 "<<endl;	

	/* no core point, but still can become a boundary point of another */
	if (seeds.size() < _minlns) 
	{
		//cout << " buildClusterSet 1 if "<<endl;	

		changeClId(Point,NOISE);		
		return false;
	
	}

	/* all points in seeds are density-reachable from Point */
	else 
	{		
		changeClIds(seeds, ClId); /* put ClId value to all members in seeds	*/
		seeds.erase(Point);	
						
		int ii = 0;
		//cout << " buildClusterSet 1 else  "<<endl;	

		while (!seeds.empty())
		{				
			ii ++;
			//cout << " buildClusterSet   WHILE ---  [ "<< ii << "] / seeds.size() "<< seeds.size()<< " / "<<  seeds.size() <<endl;
	
			int currentP = *(seeds.begin());
			set<int> result;
			//cout << " Begins region Query "<<endl;	

			regionQuery(currentP, _eps, result);
			//cout << " END region Query "<<endl;	

			
			if( result.size() >= _minlns)
			{
				for (set<int>::iterator i=result.begin(); i!=result.end(); i++)
				{
					//cout << " For buildClusterSet FOR *** i [not matches]"<<endl;	

					int resultP = *i;					
					
					if(clId(resultP)<=NOISE)
					{
						if(clId(resultP)==UNCLASSIFIED)
						{
							seeds.insert(resultP);
						}
						changeClId(resultP, ClId);						
					}
				} 
			}		
			seeds.erase(currentP);			
		} 				
		return true;
	}
}

/*---------------------------------------------------------------------
	DESCRIPTION: return seeds(objects) that are within the radius
				 from the lid point
    PARAMETERS  		
		1. lid: query point
		2. radius: query radius
    RETURN 
        1. rst: results (a set of OID)
    AUTHOR: Hoyoung Jeung, 15 June 2007    
	class PolyLine  : public vector<Point*>
*---------------------------------------------------------------------*/
inline void DbscanLine::regionQuery(int query, double r, set<int>& rst)
{    
	if(_idx && _idx->isBuilt() )// with a grid index	
		_idx->rangeQuery(query,r,rst);		
	else if(_gl && _gl->isBuilt())
		_gl->rangeQuery(query,r,rst);
	else // without index
	{			
		PolyLine* q = &_lns->at(query);  // query	
		int n=_lns->size();
		for (int i=0; i<n; i++)
		{	
			PolyLine* l = &_lns->at(i);
			//double r2 = (r + q->tol + l->tol) * (r + q->tol + l->tol);
			double r2 = r + q->tol + l->tol;
			r2 = r2*r2;
			double d2;
			q->ext.min_dist2(&l->ext,&d2);
			if( d2 <= r2)
			{
				if(_star)  // CuTS*
					q->distLL2_ST(l,&d2);
				else // CuTS
					q->distLL2(l,&d2);

				if(d2 <= r2)
					rst.insert(i);	
			}			
		}
	}	
}

/*---------------------------------------------------------------------
	DESCRIPTION: set seed points with the given id value
    PARAMETERS  
		1. seeds: seed oid set
		2. cid: setting value	
    AUTHOR: Hoyoung Jeung, 26,oct, 2005 (modified 7 June 2006)
*---------------------------------------------------------------------*/
inline void DbscanLine::changeClIds(set<int>& seeds, int cid)
{
	for (set<int>::iterator i=seeds.begin(); i!=seeds.end() ; i++)
	{		
		if(_cids[*i] <CLASSIFIED_START) //=> result is different why?
			_cids[*i]=cid;								
		//else
		//	cerr << "DbscanLine::changeClIds - overriding" << endl;
	}
}
/*------------------------------------------------------------------
	DESCRIPTION: construct the result cluster set from cid array(cids)
    AUTHOR: Hoyoung Jeung, 7 June 2006
*---------------------------------------------------------------------*/
inline void DbscanLine::buildClustSet()
{     
	int n = _lns->size();
	vector<Convoy> tmp;
	//cout << " buildClusterSet 1" <<endl;	
	for (int i=0; i<n; i++)
    {		
        int cid = _cids[i];
//		cout << " buildClusterSet 2 i["<<i<<"]"<<endl;	
		if(cid>NOISE)
        {		
			/* enlarging the storage to store new clusters */
			while(tmp.size()<cid+1)
			{
				Convoy c;
				tmp.push_back(c);
			} 
            tmp[cid].insert(_lns->at(i).o); 
        }
    }
	for (int j=0; j<tmp.size(); j++)
		if(tmp[j].size()>=_minlns)
			_clsts.push_back(tmp[j]);

	// computing each cluster's lifetime and start time
	/*			int lts = _lns->at(i).ext.mint;
	int lte = _lns->at(i).ext.maxt;

	if (_clsts[cid]._ts < 0)
	_clsts[cid]._ts = lts;
	else
	_clsts[cid]._ts = MIN(_clsts[cid]._ts,lts);

	_clsts[cid]._te = MAX(_clsts[cid]._te, lte);	*/		
}
