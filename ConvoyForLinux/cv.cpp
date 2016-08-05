#include "cv.h"

//#include "stdafx.h"
#include "cv.h"
#include "func.h"
#include <algorithm>
#include <iterator>


Convoy::Convoy(int start, int end)
{
	_ts = start;
	_te = end;
}



/*---------------------------------------------------------------------
DESCRIPTION: various operators of convoy class
PARAMETERS:
- Convoy& c : another convoy to compare with
AUTHOR: Hoyoung Jeung, 28/10/2007
*---------------------------------------------------------------------*/
bool Convoy::operator==(Convoy& c)
{
	if (c._ts != _ts || c._te != _te)
		return false;

	set<int> is;
	set_intersection(begin(), end(), c.begin(), c.end(), insert_iterator<set<int>>(is, is.begin()));
	return (is.size() == size());
}
bool Convoy::operator!=(Convoy& c)
{
	return !(*this == c);
	//bool check1 = (c._ts==_ts && c._te==_te);
	//set<int> is;
	//set_intersection(begin(),end(),c.begin(),c.end(),insert_iterator<set<int>>(is,is.begin()));	
	//bool check2 = (is.size()==size());
	//return (!check1 && !check2);
}

// check if a given convoy is a superset of this
bool Convoy::operator<=(Convoy& c)
{
	if (c._ts>_ts || c._te<_te)
		return false;

	set<int> dif;
	set_difference(begin(), end(), c.begin(), c.end(), insert_iterator<set<int>>(dif, dif.begin()));
	return (dif.size() == 0 && size() <= c.size());
}

// check if a given convoy is a subset of this
bool Convoy::operator>=(Convoy& c)
{
	if (c._ts<_ts || c._te>_te)
		return false;

	set<int> dif;
	set_difference(c.begin(), c.end(), begin(), end(), insert_iterator<set<int>>(dif, dif.begin()));
	return (dif.size() == 0 && size() >= c.size());
}

// check if a given convoy is a superset of this
bool Convoy::operator<(Convoy& c)
{
	if (c._ts >= _ts || c._te <= _te)
		return false;

	set<int> dif;
	set_difference(begin(), end(), c.begin(), c.end(), insert_iterator<set<int>>(dif, dif.begin()));
	return (dif.size() == 0 && size()<c.size());
}

// check if a given convoy is a subset of this
bool Convoy::operator>(Convoy& c)
{
	if (c._ts <= _ts || c._te >= _te)
		return false;

	set<int> dif;
	set_difference(c.begin(), c.end(), begin(), end(), insert_iterator<set<int>>(dif, dif.begin()));
	return (dif.size() == 0 && size()>c.size());
}
void Convoy::operator+=(Convoy& c)
{
	set<int> uni;
	set_union(begin(), end(), c.begin(), c.end(), insert_iterator<set<int>>(uni, uni.begin()));
	clear();
	insert(uni.begin(), uni.end());
	_ts = MIN(_ts, c._ts);
	_te = MAX(_te, c._te);
}


/*---------------------------------------------------------------------
DESCRIPTION:  check if the given oid exists during time[ts,te]
PARAMETERS
1. oid: given id to check
2. time: given time to check
RETURN
1. true: there is duplicate
2. false: no duplicate
AUTHOR: Hoyoung Jeung, 10 June 2006
*---------------------------------------------------------------------*/
bool Convoy::has(int oid, int time)
{
	// when no given data, objs.find() will return objs.end()
	bool c1 = has(oid);
	bool c2 = _ts <= time;
	bool c3 = _te >= time;
	return (c1 && c2 && c3);
}


/*---------------------------------------------------------------------
DESCRIPTION: add a member to this convoy
PARAMETERS
1. oid: object id
AUTHOR: Hoyoung Jeung, 1 June 2006
*---------------------------------------------------------------------*/
void Convoy::add(int oid)
{
	if (!insert(oid).second) // if oid already exists
		printf("Convoy::addMember - added member already exists");
}


/*---------------------------------------------------------------------
DESCRIPTION: add a member to this convoy
PARAMETERS
1. oids: a set of object id
AUTHOR: Hoyoung Jeung, 9 June 2006
*---------------------------------------------------------------------*/
void Convoy::add(set<int>& oids)
{
	for (set<int>::iterator oid = oids.begin(); oid != oids.end(); oid++)
		add(*oid);
}


/*---------------------------------------------------------------------
DESCRIPTION: show all information about this instance
AUTHOR: Hoyoung Jeung, 1 June 2006
NOTES: debugging purpose
*---------------------------------------------------------------------*/
void Convoy::print()
{
	cout << "k=" << _te - _ts + 1 << "[" << _ts << "," << _te << "], m="
		<< size() << "[";
	for (set<int>::iterator it = begin(); it != end(); it++)
	{
		if (it != begin())
			cout << ",";
		cout << *it;
	}
	cout << "]" << endl;
}



/*---------------------------------------------------------------------
DESCRIPTION: write this convoy information to a file
PARAMETERS
fp: file pointer
AUTHOR: Hoyoung Jeung, 3 June 2006
*---------------------------------------------------------------------*/
void Convoy::write(ofstream& fout)
{
	fout << "k=" << _te - _ts + 1 << "[" << _ts << "," << _te << "], m="
		<< size() << "[";
	for (set<int>::iterator it = begin(); it != end(); it++)
	{
		if (it != begin())
			fout << ",";
		fout << *it;
	}
	fout << "]" << endl;
}

/*---------------------------------------------------------------------
DESCRIPTION: show all information about a convoy set
AUTHOR: Hoyoung Jeung, 2 June 2006
NOTES: debugging purpose
*---------------------------------------------------------------------*/
void Convoy::print(vector<Convoy>& fs)
{
	//printf("number of cv: %d\n",fs->size());						
	for (int i = 0; i<fs.size(); i++)
	{
		printf("--------\n");
		fs[i].print();
	}
}
