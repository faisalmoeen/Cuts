#include "buf.h"
#include "func.h"
#include <set>
#include <fstream>
#include <iostream>
#include <limits>

Buffer::~Buffer()
{
	for (list<Point*>::iterator i = _splits.begin(); i != _splits.end(); i++)
		delete *i;
}

/*---------------------------------------------------------------------
DESCRIPTION: write the buffer to a file
PARAMETERS:
- int ts : start timestamp
- int num : the length of timestamps to write
- char* rst : file name to write
- int npoint : the number of floating points
AUTHOR: Hoyoung Jeung, 1 June 2007
*---------------------------------------------------------------------*/
void Buffer::save(int ts, int num, char* rst, int npoint)
{
	ofstream fout(rst, ios::out);
	fout.setf(ios_base::fixed, ios_base::floatfield);
	Database::iterator it = _db.begin();
	while (it != _db.end())
	{
		int oid = it->first;
		Trajectory* trj = &(it->second);
		for (int i = ts; i<ts + num; i++)
		{
			if (i >= trj->size())
				break;
			fout << oid << ' ' << trj->at(i).t << ' ' << trj->at(i).x << ' ' << trj->at(i).y << "\n";
		}
		it++;
	}
	fout.close();
}


/*---------------------------------------------------------------------
DESCRIPTION: load all data into the memory
AUTHOR: Hoyoung Jeung, 27 May 2007
NOTES: input (db file) format
OID TIME X Y
0 0 5575.52883 6390.80439
1 0 5681.62663 6054.95842
0 1 5545.73954 6356.94277
1 1 5506.18248 6167.12592
:
*---------------------------------------------------------------------*/
void Buffer::load(char *dbfile)
{
	char buf[1024];
	_numpts = 0;
	_avglen = 0;
	ifstream of(dbfile);
	if (!of)
	{
		//	cerr << "Buffer::load - invalid filename \"" << of << endl;
		return;
	}
	set<int> times;
	while (of.getline(buf, 1024) && of.good() && !of.eof())
	{
		Point p;
		int oid = p.parse(buf);
		if (oid<0)
		{
			cerr << "Buffer::load - incorrect input data at line" << _numpts << endl;
			continue;
		}
		if (_db.find(oid) == _db.end()) // no id exists
		{
			Trajectory trj;
			trj.push_back(p);
			_db.insert(Database::value_type(oid, trj));
		}
		else // already exists
			_db[oid].push_back(p);

		_ext.unionWith(p); // data extent
		times.insert(p.t);

		if (_numpts++ % 100000 == 0 && _numpts>1)
			cout << _numpts - 1 << " points have been loaded" << endl;
	}
	_ts2t.reserve(times.size());
	set<int>::iterator time = times.begin();
	int ts = 0;
	while (time != times.end())
	{
		_ts2t.push_back(*time);
		_t2ts[*time] = ts;
		time++;
		ts++;
	}
	_timedb_t = buildtimedb();

	of.close();
	summary();
}


/*---------------------------------------------------------------------
DESCRIPTION: re-arrange db according to time to speed up processing
PARAMETERS
AUTHOR: Hoyoung Jeung, 14/12/2007
*---------------------------------------------------------------------*/
clock_t Buffer::buildtimedb()
{
	clock_t  s = clock();
	long long lensum = 0;
	Database::iterator di = _db.begin();
	int maxobj = 0;
	while (di != _db.end())
	{
		Trajectory* trj = &(di->second);
		int n = trj->size();
		int oid = di->first;
		lensum += n;
		for (int i = 0; i<n; i++)
		{
			Point* p1 = &trj->at(i);
			Point* p2 = (i == n - 1) ? p1 : &trj->at(i + 1);
			if (_tlookup[p1->t].capacity()<_db.size())
				_tlookup[p1->t].reserve(_db.size());

			_tlookup[p1->t].push_back(IPoint(oid, p1));

			maxobj = MAX(_tlookup[p1->t].size(), maxobj);

			// put intermediate points between p1 and p2 if any
			int t = nextof(p1->t);
			while (t<p2->t) // split
			{
				double d = (double)(t - p1->t) / (p2->t - p1->t);
				_splits.push_back(new Point(t, p1->x + d * (p2->x - p1->x), p1->y + d * (p2->y - p1->y)));
				Point* tmp = _splits.back();
				_tlookup[t].push_back(IPoint(oid, tmp));
				t = nextof(t);

				maxobj = MAX(_tlookup[t].size(), maxobj);
			}
		}
		di++;
	}
	_avglen = (int)lensum / _db.size();
	return clock() - s;
}

void Buffer::dumpTimeDB(char* outf)
{
	ofstream fout(outf);
	fout.setf(ios_base::fixed, ios_base::floatfield); fout.precision(1);
	fout << "time : # of objects having a point at the time" << endl;
	map<int, vector<IPoint>>::iterator tb = _tlookup.begin();
	for (; tb != _tlookup.end(); tb++)
		fout << "t=" << tb->first << " : " << tb->second.size() << endl;
	fout.close();
}
/*---------------------------------------------------------------------
DESCRIPTION:
- char *inf :
- char* outf :
AUTHOR: Hoyoung Jeung, 31/10/2007
NOTE:
*---------------------------------------------------------------------*/
void Buffer::timestat(char *inf, char* outf)
{
	char buf[1024];
	ifstream fin(inf);
	map<int, int> tc;
	while (fin.getline(buf, 1024) && fin.good() && !fin.eof())
	{
		Point p;
		int oid = p.parse(buf);
		if (oid<0)
			cerr << "Buffer::load - incorrect input data\n";

		if (tc.find(p.t) == tc.end()) // no found
			tc.insert(map<int, int>::value_type(p.t, 1));
		else
			tc[p.t]++;
	}
	fin.close();

	ofstream fout(outf);
	fout.setf(ios_base::fixed, ios_base::floatfield); fout.precision(1);
	fout << "time : # of objects having a point at the time" << endl;
	for (map<int, int>::iterator i = tc.begin(); i != tc.end(); i++)
		fout << "t=" << i->first << " : " << i->second << endl;
	fout.close();
	cout << "time statistics - done" << endl;
}

/*---------------------------------------------------------------------
DESCRIPTION:
PARAMETERS:
- char *dbfile :
- char* newdbfile :
AUTHOR: Hoyoung Jeung, 27/10/2007
NOTE:
*---------------------------------------------------------------------*/
void Buffer::extractSameTime(char *dbfile, char* newdbfile)
{
	/* loading */
	char buf[1024];
	int nline = 1;
	map<int, vector<IPoint2>> tfound; // time, object count
	set<int> oids;
	ifstream of(dbfile);
	while (of.getline(buf, 1024) && of.good() && !of.eof())
	{
		vector<string> tok = Function::tokenize(buf, " \t");
		IPoint2 p;
		p.t = atoi(tok[1].c_str());
		p.x = atof(tok[2].c_str());
		p.y = atof(tok[3].c_str());
		p.oid = atoi(tok[0].c_str());

		if (oids.find(p.oid) == oids.end()) // no id exists		
			oids.insert(p.oid);

		if (tfound.find(p.t) == tfound.end()) // no found
		{
			vector<IPoint2> pts;
			pts.push_back(p);
			tfound.insert(map<int, vector<IPoint2>>::value_type(p.t, pts));
		}
		else
			tfound[p.t].push_back(p);

		if (nline++ % 10000 == 0)
			cout << nline - 1 << " points have been loaded" << endl;
	}
	of.close();
	cout << oids.size() << " objects - loaded" << endl;

	/* filtering */
	ofstream fout(newdbfile);
	fout.setf(ios_base::fixed, ios_base::floatfield);
	int found = 0;
	int n = 0;
	int nobj = oids.size();
	int ntime = tfound.size();
	map<int, vector<IPoint2>>::iterator it = tfound.begin();
	while (it != tfound.end())
	{
		int ss = it->second.size();
		if (ss == nobj)
		{
			found++;
			for (int k = 0; k<ss; k++)
			{
				IPoint2* p = &(it->second.at(k));
				fout << p->oid << " " << p->t << " " << p->x << " " << p->y << endl;
			}
		}
		it++;
		Function::progress(++n, ntime);
	}
	fout.close();
	cout << "filtering (" << tfound.size() << "->" << found << ") - done\n";
}



/*---------------------------------------------------------------------
DESCRIPTION: pre-processing of Beijing data
- char *dbfile : input
- char* newdbfile : output
AUTHOR: Hoyoung Jeung, 29/10/2007
input							output
OID(string) TIME X Y							OID(string) TIME X Y
BK142 14205327 5575.52883 6390.80439		0 0 5575.52883 6390.80439
CD102 14205331 5681.62663 6054.95842		1 0 5681.62663 6054.95842

note!! we encode the output oid as new id ignore day information to regard
differnt days' movements as those in the same day with more objects.
for example, movements of oid=BK142  => 1
*---------------------------------------------------------------------*/
void Buffer::beijing(char *dbfile, char* newdbfile)
{
	beijing(dbfile, newdbfile, -1);
}
void Buffer::beijing(char *dbfile, char* newdbfile, int numobj)
{
	char buf[1024];
	int nline = 0;
	ifstream fin(dbfile);
	map<string, Trajectory> db;
	map<string, set<int>> tfound; // oid,time found
	fin.getline(buf, 1024); // skip the first line of comments
	while (fin.getline(buf, 1024) && fin.good() && !fin.eof())
	{
		Point p;
		string oid;
		if (!p.parse(buf, oid))
		{
			cerr << "Buffer::beijing - parsing error at line " << nline << endl;
			continue;
		}
		// ex) p.t = 16035415 (2day2hour2min2sec)
		//int yy = 7; // 2007 - data does not have these information
		//int mo = 7; // July - data does not have these information	
		int dd = (int)p.t / 1000000;

		// NOTE : we consider different days movements as the same		
		int hh = (int)p.t % 1000000 / 10000;
		int mm = (int)p.t % 10000 / 100;
		int ss = (int)p.t % 100;
		//int ss = (int) p.t%100;
		//int t=ENCODE_TIME(yy,mo,dd,hh,mm,ss);int yy_,mo_,dd_,hh_,mm_,ss_;
		//yy_=YEAR(t);mo_=MONTH(t);dd_=DAY(t);hh_=HOUR(t);mm_=MINUTE(t);ss_=SECOND(t);		
		//p.t = ENCODE_TIME(0,0,dd,hh,mm,ss); // every sec
		p.t = ENCODE_TIME(0, 0, dd, hh, mm, 0);// every min
		if (db.find(oid) == db.end()) // no id exists
		{
			Trajectory trj;	trj.push_back(p);
			db.insert(map<string, Trajectory>::value_type(oid, trj));
			set<int> tt; tt.insert(p.t);
			tfound.insert(map<string, set<int>>::value_type(oid, tt));
		}
		else // already exists
		{
			if (tfound[oid].find(p.t) == tfound[oid].end()) // first
			{
				tfound[oid].insert(p.t);
				db[oid].push_back(p);
			}
		}
		if (nline++ % 100000 == 0 && nline>1)
			cout << nline - 1 << " points have been loaded" << endl;
	}
	fin.close();

	ofstream fout(newdbfile);
	fout.setf(ios_base::fixed, ios_base::floatfield); fout.precision(0);
	map<string, Trajectory>::iterator it = db.begin();
	int newoid = 1;
	int cnt = 0;
	if (numobj<0)
		numobj = db.size();
	while (it != db.end() && cnt < numobj)
	{
		Trajectory* trj = &(it->second);
		Function::sort(trj);
		for (int i = 0; i<trj->size(); i++) // time granularity => minute			
			fout << newoid << ' ' << trj->at(i).t << ' ' << trj->at(i).x << ' ' << trj->at(i).y << "\n";
		it++;
		newoid++;
		Function::progress(newoid, db.size());
		cnt++;
	}
	fout.close();
	cout << "Beijing data processing - done" << endl;
}


/*---------------------------------------------------------------------
DESCRIPTION: pre-processing of Beijing data
- char *dbfile : input
- char* newdbfile : output
AUTHOR: Hoyoung Jeung, 29/10/2007
Trucks dataset consists of 276 trajectories of 50 trucks delivering
concrete to several construction places around Athens metropolitan area
in Greece for 33 distinct days.
INPUT   (lat, lon) is in WGS84, (x, y) is in GGRS87
obj-id, traj-id, date(dd/mm/yyyy), time(hh:mm:ss), lat, lon, x, y
0862;1;10/09/2002;09:15:59;23.845089;38.018470;486253.80;4207588.10
0862;1;10/09/2002;09:16:29;23.845179;38.018069;486261.60;4207543.60
OUTPUT
OID TIME X Y
0 0 5575.52883 6390.80439
1 0 5681.62663 6054.95842
note!! we encode the output oid as obj-id+day to regard
differnt days' movements as those in the same day with more objects.
for example, movements of oid=1 at 31st => 131
*---------------------------------------------------------------------*/
void Buffer::truck(char *inf, char* outf)
{
	char buf[1024];
	ifstream fin(inf);
	ofstream fout(outf);
	fout.setf(ios_base::fixed, ios_base::floatfield); fout.precision(1);
	int lastt = 0;
	while (fin.getline(buf, 1024) && fin.good() && !fin.eof())
	{
		vector<string> tok = Function::tokenize(buf, ";");
		if (tok.size() != 8)
		{
			cout << "Buffer::truck - warning at " << buf << endl;
			continue;
		}

		vector<string> date = Function::tokenize(tok[2], "/");
		//int yy = (atoi(date[2].c_str()))%100;  // remove first two numbers like 20 of 2002
		//int mo = atoi(date[1].c_str());
		int dd = atoi(date[0].c_str());
		int oid = atoi((tok[0] + date[0]).c_str());	// <- NOTE this!	

		vector<string> time = Function::tokenize(tok[3], ":");
		int hh = atoi(time[0].c_str());
		int mm = atoi(time[1].c_str());
		int ss = atoi(time[2].c_str());
		int t = ENCODE_TIME(0, 0, 0, hh, mm, ss); // every 30 sec sampling
												  //int t=ENCODE_TIME(0,0,dd,hh,mm,0);  // every minutes sampling
												  //int yy_,mo_,dd_,hh_,mm_,ss_; 
												  //yy_=YEAR(t);mo_=MONTH(t);dd_=DAY(t);hh_=HOUR(t);mm_=MINUTE(t);ss_=SECOND(t);				
		double x = atof(tok[6].c_str());
		double y = atof(tok[7].c_str());

		if (t>lastt)  // if this 'if cluse' does not use, every 30 sec sampling
			fout << oid << ' ' << t << ' ' << x << ' ' << y << endl;

		lastt = t;
	}
	fin.close();
	fout.close();
	cout << "Truck data processing - done" << endl;
}


/*---------------------------------------------------------------------
DESCRIPTION: integrate 199 objects and extract the akta data to a file
PARAMETERS:
- char* dir : directory of the akta files
- char* outdir : output file name
- int driver : driver oid, must be in [104,321]
- int year : year to extract (must be TWO DIGITS as 99 or 01)
- int month : month to extract
- int day : day to extract
AUTHOR: Hoyoung Jeung, 28/10/2007
INPUT
DriverID,TripID,TripPositionID,X,Y,SAT,HDOP,SIGNAL,DIST,ZONE,SPEED,TIME,DATE
123,59,1,637271,6345450,9,1.3,33,100000,0,2,10:33:09,20011128
123,59,2,637271,6345450,9,1.3,33,0,0,2,10:33:10,20011128
123,59,3,637270,6345453,7,1.7,33,3,0,2,11:02:59,20011128
OUTPUT
OID TIME X Y               => see func.h for time encoding
0 0 5575.52883 6390.80439
1 0 5681.62663 6054.95842

note!! we encode the output oid as DriverID+TripID+day and ignore day
information to regard differnt days' movements as those in the same
day with more objects. for example, movements of oid=1 at 31st => 131
*---------------------------------------------------------------------*/
void Buffer::akta(char *dir, char* outf, int year, int month, int day_s, int day_e)
{
	int s_oid = 104;
	int e_oid = 320;
	int nodata[] = { 146,148,152,158,161,162,163,180,194,196,229,234,239,240,243,249,255,308 };
	set<int> shorts(nodata, nodata + 18);
	map<int, set<int>> tfound; // oid,time found

	ofstream fout(outf);
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(0);

	for (int oid = s_oid; oid <= e_oid; oid++)
	{
		if (shorts.find(oid) != shorts.end()) continue; // no driver data

														/* composing file name from 'xxxall.txt' */
		char file[300];
		sprintf(file, "%s/%dall.txt", dir, oid);
		ifstream of(file);
		if (!of)	cerr << "opening error at " << file << endl;

		char buf[1024];
		of.getline(buf, 1024); // skip the first line of comments
		while (of.getline(buf, 1024) && of.good() && !of.eof())
		{
			vector<string> tok = Function::tokenize(buf, " ,");
			if (tok.size() != 13)
			{
				cout << "Buffer::akta - warning at " << file << " : " << buf << endl;
				continue;
			}
			int date = atoi(tok[12].c_str());
			//int yy =(int) date%1000000/10000; // take only 2 numbers of a year e.g. 99=>1998, 01=>2001
			//int mo = (int) date%10000/100;
			int dd = date % 100;
			if (day_s<dd || dd>day_e)
				continue;

			vector<string> time = Function::tokenize(tok[11], ":");
			int hh = atoi(time[0].c_str());
			int mm = atoi(time[1].c_str());
			int ss = atoi(time[2].c_str());
			int t = ENCODE_TIME(0, 0, dd, hh, mm, ss);
			//int yy_,mo_,dd_,hh_,mm_,ss_;yy_=YEAR(t);mo_=MONTH(t);dd_=DAY(t);hh_=HOUR(t);mm_=MINUTE(t);ss_=SECOND(t);
			if (t>numeric_limits<long long>::max())
			{
				cerr << "Buffer::akta - file to encode time at " << file << " : " << buf << endl;
				continue;
			}
			int newid = atoi((tok[0] + tok[1]).c_str()) * 100 + dd;

			if (tfound.find(newid) == tfound.end()) // first in
			{
				set<int> tt;
				tt.insert(t);
				tfound.insert(map<int, set<int>>::value_type(newid, tt));
				fout << newid << " " << t << " " << tok[3] << " " << tok[4] << endl;
			}
			else
			{
				if (tfound[newid].find(t) == tfound[newid].end()) // first in
				{
					tfound[newid].insert(t);
					fout << newid << " " << t << " " << tok[3] << " " << tok[4] << endl;
				}
			}
		}
		of.close();
		cout << "driver " << oid << " - done" << endl;
	}

	fout.close();
	cout << "AKTA data processing - done\n";
}

// it does not integrate objects, thus will have 199 objects having longer histories
// sample every 10 seconds
void Buffer::akta2(char *dir, char* outf, int year, int month, int day_s, int day_e)
{
	int s_oid = 104;
	int e_oid = 320;
	int nodata[] = { 146,148,152,158,161,162,163,180,194,196,229,234,239,240,243,249,255,308 };
	set<int> shorts(nodata, nodata + 18);

	ofstream fout(outf);
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(0);

	for (int oid = s_oid; oid <= e_oid; oid++)
	{
		if (shorts.find(oid) != shorts.end()) continue; // no driver data

														/* composing file name from 'xxxall.txt' */
		char file[300];
		sprintf(file, "%s/%dall.txt", dir, oid);
		ifstream of(file);
		if (!of)	cerr << "opening error at " << file << endl;

		char buf[1024];
		of.getline(buf, 1024); // skip the first line of comments
		int lastt = 0;
		while (of.getline(buf, 1024) && of.good() && !of.eof())
		{
			vector<string> tok = Function::tokenize(buf, " ,");
			if (tok.size() != 13)
			{
				cout << "Buffer::akta2 - warning at " << file << " : " << buf << endl;
				continue;
			}
			int date = atoi(tok[12].c_str());
			int yy = (int)date % 1000000 / 10000; // take only 2 numbers of a year e.g. 99=>1998, 01=>2001
			int mo = (int)date % 10000 / 100;
			int dd = date % 100;

			if (yy != year) continue;
			if (mo != month) continue;
			if (dd<day_s || dd>day_e) continue;

			vector<string> time = Function::tokenize(tok[11], ":");
			int hh = atoi(time[0].c_str());
			int mm = atoi(time[1].c_str());
			int ss = atoi(time[2].c_str());
			//    ss /= 30; // 30 seconds sampling ratio
			//	ss *= 30; 
			int t = ENCODE_TIME(0, 0, dd, hh, mm, 0); // min
													  //int t=ENCODE_TIME(0,0,dd,hh,mm,ss); // sec
													  //int yy_,mo_,dd_,hh_,mm_,ss_;yy_=YEAR(t);mo_=MONTH(t);dd_=DAY(t);hh_=HOUR(t);mm_=MINUTE(t);ss_=SECOND(t);
			if (t>numeric_limits<long long>::max())
			{
				cerr << "Buffer::akta2 - file to encode time at " << file << " : " << buf << endl;
				continue;
			}

			if (t>lastt)
				fout << oid << " " << t << " " << tok[3] << " " << tok[4] << endl;

			lastt = t;
		}
		of.close();
		cout << "driver " << oid << " - done" << endl;
	}

	fout.close();
	cout << "AKTA data processing - done\n";
}

/*---------------------------------------------------------------------
DESCRIPTION: show all information about this instance
AUTHOR: Hoyoung Jeung, 28 May 2007
NOTES: debugging purpose
*---------------------------------------------------------------------*/
void Buffer::print()
{
	cout.precision(0);cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << "\n----- objects :" << _db.size() << "(e, d, t)" << endl;

	Database::iterator di = _db.begin();
	while (di != _db.end())
	{
		Trajectory* ntrj = &(di->second);

		cout << "oid=" << di->first << " :";
		for (int i = 0; i<ntrj->size(); i++)
			//cout << " (" << ntrj->at(i).e->_id << "," << ntrj->at(i).d << "," <<  ntrj->at(i).t << "),";
			cout << endl;
		di++;
	}
}

/*---------------------------------------------------------------------
DESCRIPTION: create matlab commands for visualisation
PARAMETERS:
- Trajectory& trj :
AUTHOR: Hoyoung Jeung, 29/10/2007
OUTPUT
x = [1 2 3 4 5 6 7 8 9 10];
y = [1.2 3.4 5.5 1.6 1.8 3.3 4.4 3.8 2.7 2.1];
plot(x,y)

some matlab options for drawing
y     yellow        .     point              -     solid
m     magenta       o     circle             :     dotted
c     cyan          x     x-mark             -.    dashdot
r     red           +     plus               --    dashed
g     green         *     star
b     blue          s     square
w     white         d     diamond
k     black         v     triangle (down)
*---------------------------------------------------------------------*/
void Buffer::matlab(vector<Trajectory>& trjs, char* outf, int dim, int length)
{
	ofstream fout(outf, ios::out);
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(0);

	char color[] = { 'k','r','b','g','c','m' };

	for (int k = 0; k<trjs.size(); k++)
	{
		Trajectory trj = trjs[k];
		int len = (length<0) ? trj.size() : MIN(length, trj.size());

		fout << "x" << k << "=[";
		for (int i = 0; i<len; i++)
			fout << trj[i].x << " ";
		fout << "]" << endl;

		fout << "y" << k << "=[";
		for (int i = 0; i<len; i++)
			fout << trj[i].y << " ";
		fout << "]" << endl;

		if (dim == 3)
		{
			fout << "t" << k << "=[";
			for (int i = 0; i<len; i++)
				fout << trj[i].t << " ";
			fout << "]" << endl;
		}
	}

	// ex => plot(x0,y0,'k',x1,y1,'r',x2,y2,'b'
	fout << "\nplot(";
	for (int i = 0; i<trjs.size(); i++)
	{
		int rand = i % 6; // random color
		fout << "x" << i << ",y" << i;
		if (dim == 3)
			fout << ",t" << i;
		fout << ",'" << color[i % 6] << ".-'";	 // shape of each point / line. see the option
		if (i<trjs.size() - 1)
			fout << ",";
	}
	fout << "),xlabel('X'),ylabel('Y'),";
	if (dim == 3)
		fout << "zlabel('T'),";

	// ex, legend('Trj 1','Trj 2','Trj 3','Trj 4')

	fout << "legend(";
	for (int i = 0; i<trjs.size(); i++)
	{
		fout << "'trajectory " << i << "'";
		if (i<trjs.size() - 1)
			fout << ",";
	}

	//fout << "),grid" << endl;
	fout << ")" << endl;

	fout.close();
}

// avg length of timestamps on trajectories
/*---------------------------------------------------------------------
DESCRIPTION:
PARAMETERS:
AUTHOR: Hoyoung Jeung, 1/11/2007
NOTE:
*---------------------------------------------------------------------*/
int Buffer::avglen()
{
	if (_avglen>0)
		return _avglen;
	Database::iterator di = _db.begin();
	long len = 0;
	while (di != _db.end())
	{
		len += di->second.size();
		di++;
	}
	_avglen = (int)len / _db.size();
	return _avglen;
}

/*---------------------------------------------------------------------
DESCRIPTION:
PARAMETERS:
- ofstream& fo :
AUTHOR: Hoyoung Jeung, 5/11/2007
NOTE:
*---------------------------------------------------------------------*/
void Buffer::summary(ofstream& fo)
{
	fo << "N=" << num_obj() << ", T=" << num_tstamp()
		<< ", avg trj=" << avglen() << ", num_points=" << _numpts << endl;
}

void Buffer::summary()
{
	cout << "N=" << num_obj() << ", T=" << num_tstamp()
		<< ", avg trj=" << avglen() << ", num_points=" << _numpts << endl;
}