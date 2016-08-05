#pragma once

#ifndef __CONVOY__
#define __CONVOY__

#include <set>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;


class Convoy : public set<int>
{
public:
	Convoy() { _ts = _te = -1;_assigned = false; };
	Convoy(int start, int end);
	virtual ~Convoy() {};

	int		_ts; // start of convoy duration (time)
	int		_te; // end of convoy duration (time)
	bool	_assigned;  // for moving cluster

	int		lifetime() { return _te - _ts + 1; }; // the time length of this convoy
	bool	has(int oid) { return(find(oid) != end()); };
	bool	has(int oid, int time);  // check if the given oid exists	 
	void	add(int oid);
	void	add(set<int>& oids);
	void	setLifetime(int start, int end) { _ts = start;_te = end; };

	void	print(); // display instance information to the console
	void	write(ofstream& fout);
	static void print(vector<Convoy>& sf); // information of a convoy set

	bool	operator==(Convoy& c);
	bool	operator!=(Convoy& c);
	bool	operator<=(Convoy& c); // given c is a superset of this?
	bool	operator>=(Convoy& c); // given c is a subset of this?
	bool	operator<(Convoy& c);
	bool	operator>(Convoy& c);
	void	operator+=(Convoy& c);

protected:

};

#endif 

