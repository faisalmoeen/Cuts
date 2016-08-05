/*
#pragma once
class geom
{
public:
	geom(void);
	~geom(void);
};

*/

#ifndef __GEOMETRY__
#define __GEOMETRY__

#include <vector>
#include <math.h>
#include <float.h>
#include <limits>

#ifndef MAX
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ( (x) > (y) ? (x) : (y) )
#endif
#ifndef PDIST2
#define PDIST2(p1,p2) (((p2).x-(p1).x)*((p2).x-(p1).x) + ((p2).y-(p1).y)*((p2).y-(p1).y))
#define PDIST(p1,p2)  sqrt(PDIST2((p1),(p2)))
#endif 

#ifndef DOT
#define DOT(u,v)   ((u).x * (v).x + (u).y * (v).y)
#define NORM(v)    sqrt(DOT(v,v))  // norm = length of vector
#define NORM2(v)    DOT(v,v)  // norm = length of vector
//#define SMALL_NUM  DBL_EPSILON //.2204460492503131e-016 // anything that avoids division overflow
#define SMALL_NUM  0.0000001 //.2204460492503131e-016 // anything that avoids division overflow
#define LARGE_NUM  DBL_MAX //.2204460492503131e-016 // anything that avoids division overflow
#endif 

#ifndef PI
#define PI 3.14159265358979323846
#endif 
#ifndef ABS
#define ABS(i) 	( (i) >= 0 ? (i) : -(i) )
#endif

#define WRONG_VALUE (PI)

using namespace std;


class Point
{
public:
	Point(){t=-1;}
	Point(char* buf){parse(buf);}
	Point(double xx, double yy){x=xx;y=yy;}
	Point(int tt, double xx, double yy){t=tt;x=xx;y=yy;}
	virtual ~Point(){};

	int		t;
	double	x;
	double	y;

	bool operator ==(const Point& p){return(t==p.t&&x==p.x&&y==p.y);}
	bool operator !=(const Point& p){return !(*this==p);}
	Point Point::operator+( const Point& p){return Point(x+p.x,y+p.y);}
	Point Point::operator+( const double c){return Point(x+c,y+c);}
	Point Point::operator-( const Point& p ){return Point(x-p.x,y-p.y);}
	Point Point::operator*( const double c ){return Point(x*c,y*c);} 
	Point Point::operator/( const double c ){return Point(x/c,y/c);} 
	int  parse(char* buf);
	bool parse(char* buf, string& oid);
};

class Point2D
{
public:
	Point2D(){};
	Point2D(double xx, double yy){x=xx;y=yy;}
	virtual ~Point2D(){};

	double	x;
	double	y;

	Point2D Point2D::operator=( const Point& p){return Point2D(x=p.x,y=p.y);}
	Point2D Point2D::operator+( const Point2D& p){return Point2D(x+p.x,y+p.y);}
	Point2D Point2D::operator-( const Point2D& p ){return Point2D(x-p.x,y-p.y);}
	Point2D Point2D::operator+( const Point& p){return Point2D(x+p.x,y+p.y);}
	Point2D Point2D::operator-( const Point& p ){return Point2D(x-p.x,y-p.y);}
	Point2D Point2D::operator*( const double c ){return Point2D(x*c,y*c);} 
	Point2D Point2D::operator/( const double c ){return Point2D(x/c,y/c);} 
};

struct IPoint
{
	int		o; // object id
	Point*	p; // point pointer
	IPoint(int oid, Point * pt){o=oid;p=pt;};
	IPoint(){};
};

class PolyLine;
class MBR
{
public:
	MBR(){haspoint = false;};
	MBR(double smallx, double smally, double bigx, double bigy)
	{ minx=smallx;miny=smally;maxx=bigx;maxy=bigy;haspoint=true;}
	MBR(Point2D& min, Point2D& max)
	{ minx=min.x;miny=min.y;maxx=max.x;maxy=max.y;haspoint=true;}
	MBR(Point& min, Point& max)
	{ minx=min.x;miny=min.y;maxx=max.x;maxy=max.y;haspoint=true;}
	virtual ~MBR(){};
	
	double	minx;
	double	miny;
	double	maxx;
	double	maxy;
	
	int		size(){return 4;};
	void	unionWith(Point& pt);
	void	unionWith(MBR& mbr);
	void	unionWith(Point* pt){unionWith(*pt);};
	void	unionWith(vector<Point>& pts);
	bool	contain(Point *pt);
	double	area(){	return (maxx-minx) * (maxy-miny);};
	double  width(){return maxx-minx;};
	double  height(){return maxy-miny;};
	double  diagonal(){return sqrt((maxy-miny)*(maxy-miny)+(maxx-minx)*(maxx-minx));}
	Point	center(){return Point(-1,(minx+maxx)/2,(miny+maxy)/2);}

	Point	min() {return Point(minx,miny);}
	Point	max() {return Point(maxx,maxy);}

	bool	intersect( Point &p1, Point& p2 )
	{
		double st,et,fst = 0,fet = 1;   
		double const *bmin = &minx;   
		double const *bmax = &maxx;   
		double const *si = &p1.x;  		
		double const *ei = &p2.x;

		for (int i = 0; i < 2; i++) 
		{      
			if (*si < *ei) 
			{         
				if (*si > *bmax || *ei < *bmin)
					return false;
				double di = *ei - *si;
				st = (*si < *bmin)? (*bmin - *si) / di: 0;
				et = (*ei > *bmax)? (*bmax - *si) / di: 1;
			}      
			else {
				if (*ei > *bmax || *si < *bmin)
					return false;         
				double di = *ei - *si;         
				st = (*si > *bmax)? (*bmax - *si) / di: 0;         
				et = (*ei < *bmin)? (*bmin - *si) / di: 1;      
			}      
			if (st > fst) fst = st;
			if (et < fet) fet = et;
			if (fet < fst) return false;
			bmin++; bmax++; si++; ei++;
		}  
		return true;
	}

	void MBR::min_dist( MBR &m, double* rst ) {min_dist(&m,rst);}
	void MBR::min_dist( MBR *m, double* rst ) {double d;min_dist2(m,&d);*rst=sqrt(d);}
	void MBR::min_dist2( MBR &m, double* rst ) {min_dist2(&m,rst);}
	void MBR::min_dist2( MBR* m, double* rst )
	{
		double xd,yd;
		if(maxx < m->minx)
			xd = m->minx - maxx;
		else if(minx > m->maxx)
			xd = minx - m->maxx;
		else // overlapping
			xd = 0;
		if(maxy < m->miny)
			yd = m->miny - maxy;
		else if(miny > m->maxy)
			yd = miny - m->maxy;
		else // overlapping
			yd = 0;
		*rst = xd*xd + yd*yd;				
	}

	void MBR::max_dist2( MBR* m, double* rst )
	{
		double xd = MAX(maxx,m->maxx) - MIN(minx,m->minx);
		double yd = MAX(maxy,m->maxy) - MIN(miny,m->miny);			
		*rst = xd*xd + yd*yd;				
	}

protected:
	bool	haspoint; // if this MBR does not contain any point, it is false
};


class MBR3 : public MBR
{
public:
	MBR3(){haspoint = false;};
	virtual ~MBR3(){};

	int		mint;	
	int		maxt;
	void	unionWith(Point &p)
	{
		if(haspoint)
		{	
			mint = MIN(mint, p.t);
			maxt = MAX(maxt, p.t);
			minx = MIN(minx, p.x);
			maxx = MAX(maxx, p.x);
			miny = MIN(miny, p.y);
			maxy = MAX(maxy, p.y);
		}
		else
		{
			mint = maxt	= p.t;
			minx = maxx = p.x;
			miny = maxy = p.y;
			haspoint = true;
		}
	}
	void unionWith(MBR3& mbr)
	{
		if(haspoint)
		{	
			minx = MIN(minx,mbr.minx);
			miny = MIN(miny,mbr.miny);
			mint = MIN(mint,mbr.mint);
			maxx = MAX(maxx,mbr.maxx);
			maxy = MAX(maxy,mbr.maxy);
			maxt = MAX(maxt,mbr.maxt);
		}
		else
		{
			mint = mbr.mint;
			minx = mbr.minx;
			miny = mbr.miny;
			maxt = mbr.maxt;
			maxx = mbr.maxx;
			maxy = mbr.maxy;
			haspoint = true;
		}
	}
	int		timewidth(){return maxt-mint;}
	Point	min() {return Point(mint,minx,miny);}
	Point	max() {return Point(maxt,maxx,maxy);}

private:
	bool	haspoint; // if this MBR does not contain any point, it is false
};


class PolyLine  : public vector<Point*>
{ 
public:
	int		o; 
	double  tol;
	MBR		ext;

	PolyLine(){};
	PolyLine(int oid){o=oid;};
	void	push_back(Point* p){ext.unionWith(p);vector<Point*>::push_back(p);}
	void	push_back(int oid, Point* p1, Point* p2){o=oid;push_back(p1);push_back(p2);}
	MBR*	mbr(){return &ext;}
	void	computeExtent(){for(int i=0; i<size();i++)ext.unionWith(at(i));}
	double	totalLength(){double len=0;	if(size()>1) 
								for(int i=0; i<size()-1;i++)
									len+=PDIST(*(at(i)),(*(at(i+1))));
								return len;}
	double	avgLength(){return totalLength()/size();}
	double  maxSegLength(){double maxlen=0;	if(size()>1) 
							for(int i=0; i<size()-1;i++)
								maxlen = MAX(maxlen,PDIST(*(at(i)),(*(at(i+1)))));
							return maxlen;}
	int		startTime(){return at(0)->t;}
	int		endTime(){return at(size()-1)->t;}

	static Point*	split(int t, Point* p1, Point* p2){	
		double d = (p2->t==p1->t)? 0 : (double)(t-p1->t)/(p2->t - p1->t);          
		return new Point(t,p1->x + d * (p2->x - p1->x),p1->y + d * (p2->y - p1->y));
	}
	// distance between polylines
	void distLL( PolyLine &l, double* rst ) {distLL(&l, rst);}
	void distLL( PolyLine *l, double* rst ) {double d;distLL2(l,&d); *rst = sqrt(d);}
	void distLL2( PolyLine &l, double* rst ) {distLL2(&l,rst);}
	void distLL2( PolyLine* l, double* rst )
	{
		*rst=LARGE_NUM;
		double d;
		int n1 = size()-1;
		int n2 = l->size()-1;
		int lastj=0;
		for(int i=0; i<n1; i++)
			for(int j=lastj; j<n2; j++)
			{
				if(l->at(j+1)->t <at(i)->t)	break;

				if(at(i+1)->t < l->at(j)->t) 
					continue;
				else
					lastj=j;	
				dist_ll2(at(i),at(i+1),l->at(j),l->at(j+1),&d);
				*rst = MIN(*rst,d);	 
			}
	}
	
	// distance between line segments
	void dist_ll2(Point& p1, Point& p2, Point& q1, Point& q2, double* rst) {dist_ll2(&p1,&p2,&q1,&q2,rst);}
	void dist_ll2(Point* p1, Point* p2, Point* q1, Point* q2, double* rst)
	{
		Point2D u,v,w, dP;
		u = *p2 - *p1;
		v = *q2 - *q1;
		w = *p1 - *q1; 
		double    a = DOT(u,u);        // always >= 0
		double    b = DOT(u,v);
		double    c = DOT(v,v);        // always >= 0
		double    d = DOT(u,w);
		double    e = DOT(v,w);
		double    D = a*c - b*b;       // always >= 0
		double    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
		double    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

		// compute the line parameters of the two closest points
		if (D < SMALL_NUM) { // the lines are almost parallel
			sN = 0.0;        // force using point P0 on segment S1
			sD = 1.0;        // to prevent possible division by 0.0 later
			tN = e;
			tD = c;
		}
		else {                // get the closest points on the infinite lines
			sN = (b*e - c*d);
			tN = (a*e - b*d);
			if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
				sN = 0.0;
				tN = e;
				tD = c;
			}
			else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
				sN = sD;
				tN = e + b;
				tD = c;
			}
		}

		if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
			tN = 0.0;
			// recompute sc for this edge
			if (-d < 0.0)
				sN = 0.0;
			else if (-d > a)
				sN = sD;
			else {
				sN = -d;
				sD = a;
			}
		}
		else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
			tN = tD;
			// recompute sc for this edge
			if ((-d + b) < 0.0)
				sN = 0;
			else if ((-d + b) > a)
				sN = sD;
			else {
				sN = (-d + b);
				sD = a;
			}
		}
		// finally do the division to get sc and tc
		sc = (ABS(sN) < SMALL_NUM ? 0.0 : sN / sD);
		tc = (ABS(tN) < SMALL_NUM ? 0.0 : tN / tD);

		// get the difference of the two closest points
		dP = w + (u * sc) - (v * tc);  // = S1(sc) - S2(tc)
		*rst = NORM2(dP);
		if(*rst>0.0 && *rst<SMALL_NUM)
			*rst = 0.0;
	}

	void distLL_ST( PolyLine &l, double* rst) {distLL_ST(&l,rst);}
	void distLL_ST( PolyLine *l, double* rst ) {double d;distLL2_ST(l,&d);*rst=sqrt(d);}
	void distLL2_ST( PolyLine &l, double* rst ) {distLL2_ST(&l,rst);}
	void distLL2_ST( PolyLine* l, double* rst )
	{
		double d;
		*rst=LARGE_NUM;
		int n1 = size()-1;
		int n2 = l->size()-1;
		int lastj=0;
		for(int i=0; i<n1; i++)
			for(int j=lastj; j<n2; j++)
			{
				if(l->at(j+1)->t <at(i)->t)
					break;

				if(at(i+1)->t < l->at(j)->t)
					continue;
				else
					lastj=j;

				cpa_distance2(at(i),at(i+1),l->at(j),l->at(j+1),&d);
				*rst = MIN(*rst,d);
			}
	}

	// distance of closest point of approach (CPA)
	void 	cpa_distance2( Point* p1, Point* p2, Point* q1, Point* q2, double* rst)
	{
		double  ct; // cpa time
		int		dt1 = (p2->t-p1->t);
		int		dt2 = (q2->t-q1->t);
		Point2D v1(.0,.0),v2(.0,.0),dv;
		if(dt1>0)
			v1 =  (*p2 - *p1)/(double)dt1; // velocity of p1-p2
		if(dt2>0)
			v2 = (*q2 - *q1)/(double)dt2; // velocity of q1-q2
		dv = v1 - v2;
		double  dv2 = DOT(dv,dv);

		int mint = MAX(p1->t,q1->t); // minimum overlapping time
		int maxt = MIN(p2->t,q2->t); // maximum overlapping time

		if (dv2 < SMALL_NUM)      // the tracks are almost parallel, not to have an overflow
			ct= mint;            // any time is ok.  Use time 0.
		else
		{
			Point2D  w0; w0 = *p1 - *q1;
			ct = (-DOT(w0,dv)) / (double)dv2;
		}

		if (ct<mint)
			ct = mint;
		else if (ct>maxt)
			ct = maxt;

		*rst=PDIST2((v1 * (ct-p1->t)+ *p1),((v2 * (ct-q1->t)+*q1)));           // distance at CPA
		if(*rst>0 && *rst<SMALL_NUM)
			*rst = 0.0;
	}
};

#endif


