/*
#include "stdafx.h"
#include "dbscan.h"


dbscan::dbscan(void)
{
}


dbscan::~dbscan(void)
{
}
*/

#include "stdafx.h"
#include "dbscan.h"
#include "func.h"
#include <iostream>

Dbscan::Dbscan(vector<IPoint>* ptset, int minpts, double eps)
{ 
	_pts=ptset;
	_minpts=minpts; 
	_eps=eps;
	_grid=NULL;
}

/*---------------------------------------------------------------------
	DESCRIPTION: perform dbscan algorithm
    AUTHOR: Hoyoung Jeung, 7 June 2006
*---------------------------------------------------------------------*/
void Dbscan::perform()
{
	/* initializing clusterID values */
	_cids.resize(_pts->size(),UNCLASSIFIED);   

	dbscan();
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
void Dbscan::dbscan()
{		
	int ClusterId = nextId(NOISE);	// the first cluster id begins from 0

	for (int i=0; i<_pts->size(); i++)
	{
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
inline bool Dbscan::expandCluster(int Point, int ClId)
{
	set<int> seeds; 
	regionQuery(Point, _eps, seeds);		

	/* no core point */
	if (seeds.size() < _minpts) 
	{
		changeClId(Point,NOISE);		
		return false;
	}

	/* for core points */
	else 
	{		
		changeClIds(seeds, ClId); /* put ClId value to all members in seeds	*/
		seeds.erase(Point);	
				
		while (!seeds.empty())
		{				
			int currentP = *(seeds.begin());
			set<int> result;
			regionQuery(currentP, _eps, result);
			
			if( result.size() >= _minpts)
			{
				for (set<int>::iterator i=result.begin(); i!=result.end(); i++)
				{
					int resultP = *i;					
					
					if((clId(resultP)==UNCLASSIFIED)||(clId(resultP)==NOISE))
					{
						if(clId(resultP)==UNCLASSIFIED)
						{
							seeds.insert(resultP);
						}
						changeClId(resultP, ClId);						
					}
				} // for
			} // if			
			seeds.erase(currentP);			
		} // while		
		
		return true;
	}
}

/*---------------------------------------------------------------------
	DESCRIPTION: return seeds(objects) that are within the radius
				 from the center point
    PARAMETERS  		
		1. center: query point
		2. radius: query radius
		3. rst: results (a set of OID)
    AUTHOR: Hoyoung Jeung, 26,oct, 2005 (modified Oct 2007)    
*---------------------------------------------------------------------*/
inline void Dbscan::regionQuery(int query, double radius, set<int>& rst)
{		
	if(_grid && _grid->isBuilt())// with a grid index	
		_grid->rangeQuery(_pts,query,radius,rst);		
	else // without index
	{			
		double radius2 = radius*radius;
		for (int i=0; i<_pts->size(); i++)
			if(PDIST2(*(_pts->at(query).p),*(_pts->at(i).p))<=radius2)
				rst.insert(i);
	}	
}


/*---------------------------------------------------------------------
	DESCRIPTION: set seed points with the given id value
    PARAMETERS  
		1. seeds: seed oid set
		2. cid: setting value	
    AUTHOR: Hoyoung Jeung, 26,oct, 2005 (modified 7 June 2006)
*---------------------------------------------------------------------*/
inline void Dbscan::changeClIds(set<int>& seeds, int cid)
{
	for (set<int>::iterator i=seeds.begin(); i!=seeds.end() ; i++)
	{					
		if (_cids[*i] < CLASSIFIED_START ) 
			_cids[*i]=cid;		
	}
}

/*---------------------------------------------------------------------
	DESCRIPTION: construct the result cluster set from cid array(cids)
    AUTHOR: Hoyoung Jeung, 7 June 2006
*---------------------------------------------------------------------*/
void Dbscan::buildClustSet()
{	
	//verify(); // debugging purpose
	vector<Convoy> tmp;
	for (int pid=0; pid<_pts->size(); pid++)
	{		
		int cid = _cids[pid];
		if(cid>NOISE)
		{		
			/* enlarging the space to store new clusters */
			while(tmp.size()<cid+1)
			{
				Convoy c;
				tmp.push_back(c);
			}

			int oid = _pts->at(pid).o;
			tmp.at(cid).insert(oid);			
		}
	}

	for (int j=0; j<tmp.size(); j++)
		if(tmp[j].size()>=_minpts)
			_clsts.push_back(tmp[j]);
}

