/*
#pragma once
class func
{
public:
	func(void);
	~func(void);
};

*/
#ifndef __HYFUNCTION__
#define __HYFUNCTION__

#include <CSTDLIB>
#include "geom.h"
#include <time.h>
#include <vector>
#include <map>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#ifndef PI
#define PI 3.14159265358979323846
#endif 

#ifndef M_PI
#define M_PI PI
#endif 

#ifndef PDIST2
#define PDIST2(p1,p2) (((p2).x-(p1).x)*((p2).x-(p1).x) + ((p2).y-(p1).y)*((p2).y-(p1).y))
#define PDIST(p1,p2)  sqrt(PDIST2((p1),(p2)))
#endif 

// dot product (3D) which allows vector operations in arguments
#define DOT(u,v)   ((u).x * (v).x + (u).y * (v).y)
#define NORM(v)    sqrt(DOT(v,v))  // norm = length of vector
#define SMALL_NUM  0.00000001 // anything that avoids division overflow


#ifndef MAX
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ( (x) > (y) ? (x) : (y) )
#endif

#ifndef ABS
#define ABS(i) 	( (i) >= 0 ? (i) : -(i) )
#endif

#ifndef YEAR
// time encoding 1 : change to seconds
#define YEAR(t_)	((int)((t_)/31104000))			//3600*24*30*12
#define MONTH(t_)	((int)(((t_)%31104000)/2592000))	//3600*24*30
#define DAY(t_)		((int)(((t_)%2592000)/86400))		//3600*24*30/3600*24
#define HOUR(t_)	((int)(((t_)%86400)/3600))	
#define MINUTE(t_)	((int)(((t_)%3600)/60))
#define SECOND(t_)	((int)((t_)%60))
#define ENCODE_TIME(y_,mo_,d_,h_,m_,s_)	((y_)*31104000+(mo_)*2592000+(d_)*86400+(h_)*3600+(m_)*60+(s_))

// time encoding 2 : 1116075935 -> mon(11) + day(16) + hour(07) + min(59) + sec(35)
#define MONTH1(t_)	((int)((t_)/100000000))
#define DAY1(t_)	((int)((t_)%100000000/1000000))
#define HOUR1(t_)	((int)((t_)%1000000/10000))
#define MINUTE1(t_)	((int)((t_)%10000/100))
#define SECOND1(t_)	((int)((t_)%100))
#define ENCODE_TIME1(mo_,d_,h_,m_,s_)	((mo_)100000000+(d_)*1000000+(h_)*10000+(m_)*100+(s_))
#endif

class LineSegment;
class Point;


class Function  
{
public:
	Function();
	virtual ~Function();

	/*
	 *	random functions
	 */
	static double	drand48();
	static double	uniform(double _min, double _max);
	static int		uniform(int _min, int _max);
	static int		uniform(int _min, int _max, int avoid);
	static double	new_uniform(int _d_num);
	static double	gaussian (double mean, double sigma);
	static double	zipf(double x1, double x2, double p);
	static double	new_uniform(double _min, double _max);

	/*	Moving Point functions */
	static double	velocity(Point *pts, int reftm, int retro);

	/* Misc */
	static void		progress(int now, int total);
	static int		numLine(char *fname);
	void			extent(char *fn, int nline,
					  double *maxx, double *maxy, double *minx, double *miny);
	static void		sort(vector<Point>& trj){sort(&trj);};
	static void		sort(vector<Point>* trj);

	/* string control */
	static void trim(string& str){
		string::size_type pos = str.find_last_not_of(' ');
		if(pos != string::npos) {
			str.erase(pos + 1);
			pos = str.find_first_not_of(' ');
			if(pos != string::npos) str.erase(0, pos);
		}
		else str.erase(str.begin(), str.end());};

	static bool isSame(string& s1, string& s2)
		{trim(s1);trim(s2);return (s1.compare(s2)==0);};

	static bool isSame(char* s1, char* s2)
		{return isSame(string(s1),string(s2));};

	static vector<string> tokenize(string& str, char* delim){
		vector<string> tokens;
		string deli(delim);
		// Skip delimiters at beginning.
		string::size_type lastPos = str.find_first_not_of(deli, 0);
		// Find first "non-delimiter".
		string::size_type pos     = str.find_first_of(deli, lastPos);
		while (string::npos != pos || string::npos != lastPos){
			// Found a token, add it to the vector.
			tokens.push_back(str.substr(lastPos, pos - lastPos));
			// Skip delimiters.  Note the "not_of"
			lastPos = str.find_first_not_of(deli, pos);
			// Find next "non-delimiter"
			pos = str.find_first_of(deli, lastPos);
		}
		return tokens;};

	static vector<string> tokenize(char* s, char* delim)	
		{return tokenize(string(s),delim);};
	static const char* combine(char* s1, char* s2){
		return (string(s1) + string(s2)).c_str();}

};

#endif