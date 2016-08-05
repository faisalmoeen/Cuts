/*
#include "stdafx.h"
#include "geom.h"


geom::geom(void)
{
}


geom::~geom(void)
{
}
*/
#include "stdafx.h"
#include "geom.h"
#include "func.h"
#include <fstream>
#include <iostream>



/*---------------------------------------------------------------------
	DESCRIPTION: parse given string and construct a point
	PARAMETERS: 
		- char* buf : buffer
	AUTHOR: Hoyoung Jeung, 27/10/2007 
	NOTES: input (db file) format
	OID TIME X Y
	0 0 5575.52883 6390.80439
	1 0 5681.62663 6054.95842
	0 1 5545.73954 6356.94277
	1 1 5506.18248 6167.12592
*---------------------------------------------------------------------*/
int Point::parse(char* buf)
{
	vector<string> tok = Function::tokenize(buf," \t");
	if(tok.size()!=4)
		return -1;

	t = atoi(tok[1].c_str());
	x = atof(tok[2].c_str()); 
	y = atof(tok[3].c_str());		
	return  atoi(tok[0].c_str()); 
}

bool Point::parse(char* buf, string& oid)
{
	vector<string> tok = Function::tokenize(buf," \t");
	if(tok.size()!=4)
		return false;

	t = atoi(tok[1].c_str());
	x = atof(tok[2].c_str()); 
	y = atof(tok[3].c_str());		
	oid = tok[0];
	return true;; 
}


/*---------------------------------------------------------------------
	DESCRIPTION: unionWith current mbr with the given point
    PARAMETERS  
		1. pt: a moving point to be unioned		
    RETURN 
        1. void
    AUTHOR: Hoyoung Jeung, 7,nov, 2005
    NOTES: 
*---------------------------------------------------------------------*/
inline void MBR::unionWith(Point &p)
{
	if(haspoint)
	{	
		minx = MIN(minx, p.x);
		maxx = MAX(maxx, p.x);
		miny = MIN(miny, p.y);
		maxy = MAX(maxy, p.y);
	}
	else
	{
		minx = maxx = p.x;
		miny = maxy = p.y;
		haspoint = true;
	}
}

/*---------------------------------------------------------------------
	DESCRIPTION: unionWith current mbr with the given points
    PARAMETERS  
		1. pts: moving points to be unioned
		2. nPts: the number of the given points
    RETURN 
        1. void
    AUTHOR: Hoyoung Jeung, 7,nov, 2005
    NOTES: 
*---------------------------------------------------------------------*/
void MBR::unionWith(vector<Point>& pts)
{
	for (int i=0; i< pts.size(); i++)
		unionWith(pts[i]);
}

void MBR::unionWith(MBR& mbr)
{
	if(haspoint)
	{	
		minx = MIN(minx,mbr.minx);
		miny = MIN(miny,mbr.miny);
		maxx = MAX(maxx,mbr.maxx);
		maxy = MAX(maxy,mbr.maxy);
	}
	else
	{
		minx = mbr.minx;
		miny = mbr.miny;
		maxx = mbr.maxx;
		maxy = mbr.maxy;
		haspoint = true;
	}
}

/*---------------------------------------------------------------------
	DESCRIPTION: check if the given point is inside of this mbr
    PARAMETERS  
		1. pt: a point to check
    RETURN 
        1. true: the point is inside of the mbr boundary
		2. false: outside
    AUTHOR: Hoyoung Jeung, 7,nov, 2005
    NOTES: 
*---------------------------------------------------------------------*/
inline bool MBR::contain(Point *pt)
{
	return ((pt->x>=minx) && (pt->x<=maxx) && (pt->y>=miny) && (pt->y<=maxy));
}
