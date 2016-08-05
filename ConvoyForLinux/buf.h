#pragma once

#ifndef __BUFFER__
#define __BUFFER__

#include <stdio.h>
#include <stdlib.h>
#include "geom.h"
#include <vector>
#include <set>
#include <map>
#include <time.h>
#include <list>

using namespace std;

//typedef vector<Point> Trajectory;
//typedef map<int,Trajectory> Database; // <oid, trajectory>

class Trajectory : public vector<Point>
{
public:
	void	print();
		bool operator !=(const Trajectory& t){return !(*(this)==t);}
	bool operator ==(const Trajectory& t){
		int n = size();
		if(n!=t.size()) return false;
		for(int i=0; i<n;i++)
			if(at(i)!=t[i])
				return false;
		return true;
	}

};


class Database : public map<int,Trajectory>
{
public:
	Trajectory first(){ return begin()->second;} // first trajectory
};


class Buffer 
{
public:	
	Buffer(){};
	Buffer(char *dbfile){open(dbfile);}
	virtual ~Buffer();

	Database* getDB() {return &_db;};

	void	open(char *dbfile){_avglen=0;load(dbfile);}

	MBR		extent(){return _ext;};

	long long	num_points(){return _numpts;}
	/* number of objects in DB */
	int		num_obj() {return _db.size();};

	/* maximum number of timestamps in DB */
	int		num_tstamp() {return _tlookup.size();}

	int		num_time() {return _ext.timewidth();}

	int		timeof(int timestamp){return _ts2t[timestamp];}
	int		tstampof(int time){return _t2ts[time];}

	int		nextof(int time){		
				map<int,int>::iterator ct= _t2ts.upper_bound(time);
				return (time==lasttime()) ? lasttime() : ct->first;}
	int		beforeof(int time){
				map<int,int>::iterator ct= _t2ts.lower_bound(time);
				return (time==begintime()) ? begintime() :ct->first;}

	int		lasttime() {return _ext.maxt;};

	int     begintime(){return _ext.mint;}

	int		avglen();

	void	summary(); 
	void	summary(ofstream& fo); 

	map<int,vector<IPoint>>* tlookup(){return &_tlookup;}
	void dumpTimeDB(char* outf);

	/* read object files as the RamBuffer size */
	vector<IPoint>*	objAtTstamp(int tstamp){return &_tlookup[_ts2t[tstamp]];}
	vector<IPoint>*	objAtTime(int time){return &_tlookup[time];}
	
	void	print();
	void    save(int ts, int e, char* fout, int npoint); // n point is # of points after decimal point

	static void extractSameTime(char *dbfile, char* outf);
	static void	beijing(char *dbfile, char* outf);
	static void beijing(char *dbfile, char* newdbfile, int numobj);
	static void	truck(char *dbfile, char* outf);	
	static void akta(char *dir, char* outf, int year, int month, int day_s, int day_e);	
	static void akta2(char *dir, char* outf, int year, int month, int day_s, int day_e);

	static void timestat(char *inf, char* outf);
	static void matlab( vector<Trajectory>& trjs, char* outf){matlab(trjs,outf,2,-1);};
	static void matlab( vector<Trajectory>& trjs, char* outf, int dim, int length);

	class IPoint2 : public Point
	{ public:
		int		oid;
	};

	clock_t	 _timedb_t;

protected:
	
	Database _db;
	MBR3	 _ext;	// extent of the map
	int		 _avglen;
	long long _numpts;

	map<int,vector<IPoint>> _tlookup; //<time, (oid,point pointer)>
	vector<int> _ts2t; // timestamp to time
	map<int,int> _t2ts; // time to timestamp
	/* load all data into the memory */
	void    load(char *dbfile);
	clock_t buildtimedb();
	list<Point *> _splits;

	
};

#endif

class buf
{
public:
	buf();
	~buf();
};
