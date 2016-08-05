/*
#include "stdafx.h"
#include "cvq.h"


cvq::cvq(void)
{
}


cvq::~cvq(void)
{
}
*/

#include "stdafx.h"
#include "cvq.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <cassert>
#include "func.h"
#include <iterator>

ConvoyQuery::ConvoyQuery()
{
	_buf=NULL;
	_showprogress = true;
	_e = -1.0;
	_index = false; 
	_index_t=0;
	initialize();
}


/*---------------------------------------------------------------------
	DESCRIPTION: report params/result cv information to the given file
    PARAMETERS  
		1. fname: report file name to write    
    AUTHOR: Hoyoung Jeung, 6 June 2006
    NOTES: debugging purpose
*---------------------------------------------------------------------*/
void ConvoyQuery::write(char* fname)
{
	ofstream fout(fname, ios::out); 
	fout.setf(ios_base::boolalpha);
	fout.precision(1);
	fout.setf(ios_base::fixed,ios_base::floatfield);	
	write(fout);
	fout.close();
}


/*---------------------------------------------------------------------
   DESCRIPTION: initialize all the variable to perform new query processing
   AUTHOR: Hoyoung Jeung, 25/10/2007
*---------------------------------------------------------------------*/
void ConvoyQuery::initialize()
{
	_rst.clear();	
	_elapsed_t=0;
	_fneg=_fpos=_exac=_corr=0;	
}


/*---------------------------------------------------------------------
   DESCRIPTION: 
   PARAMETERS
	  - ofstream& fo :
   AUTHOR: Hoyoung Jeung, 14/12/2007
*---------------------------------------------------------------------*/
void ConvoyQuery::write( ofstream& fo )
{
	summary(fo);
	//for(int i=0; i<_rst.size(); i++)
	//	_rst[i].write(fo);

	for(int i=0; i<_rst.size(); i++)
	{
		fo << "k=" << _rst[i]._te-_rst[i]._ts+1 << "[" 
			<< _rst[i]._ts << "," <<_rst[i]._te <<"](["
			<< _buf->tstampof(_rst[i]._ts) << "," <<_buf->tstampof(_rst[i]._te) <<"]), m="
			<< _rst[i].size() << "[";		
		for(set<int>::iterator it = _rst[i].begin(); it != _rst[i].end(); it++)
		{	
			if(it != _rst[i].begin())
				fo << ",";
			fo << *it ;
		}	
		fo << "]" << endl;
	}	
}	


/*---------------------------------------------------------------------
    DESCRIPTION: compute the number of false negatives and false positives
    PARAMETERS  
        1. ans : a convoy set which is supposed to be
    AUTHOR: Hoyoung Jeung, 20 June 2007
*---------------------------------------------------------------------*/
void ConvoyQuery::validate( vector<Convoy>& ans)
{
	_fpos = _fneg = _exac = _corr = 0;
	if(_rst.size()==0)
	{
		_fneg = ans.size();
		return;
	}
	if(ans.size()==0)
	{
		_fpos=_rst.size();
		return;
	}	
	for(vector<Convoy>::iterator a=ans.begin(); a!=ans.end(); a++)
	{
		for(vector<Convoy>::iterator c=_rst.begin(); c!=_rst.end(); c++)
		{
			//cout << "\nanswer "; a->print(); cout << "candidate "; c->print();	
			if(*a==*c)
				_exac++;
			else if(*a<=*c)	
				_corr++;
		}   
	}
	_fneg2 = (1.0 - (_exac+_corr)/(double)ans.size())*100;
	_fpos2 = (1.0 - (_exac+_corr)/(double)_rst.size())*100;
	_fneg = (1.0 - (_exac)/(double)ans.size())*100;
	_fpos = (1.0 - (_exac)/(double)_rst.size())*100;
}

/*---------------------------------------------------------------------
	DESCRIPTION: compute the range Eps for DBSCAN
	PARAMETERS: 
		- int m : minimum number of objects
	RETURNS:
        - double
	AUTHOR: Hoyoung Jeung, 22/11/2007 
	NOTE: 
*---------------------------------------------------------------------*/
double ConvoyQuery::findEps( int m, int k )
{	
	int ntrial = 300;
	ntrial = (ntrial>_buf->num_tstamp()) ? _buf->num_tstamp() : ntrial;
	int stop = 1000;
	int nsample = 500;

	set<double> es;

	map<int,vector<IPoint>>* tg = _buf->tlookup();
	multimap<int,int> tt; // (# of points, time) -> this is sorted by the number
	for(map<int,vector<IPoint>>::iterator it=tg->begin(); it!=tg->end(); it++)
		tt.insert(multimap<int,int>::value_type(it->second.size(),it->first));

	multimap<int,int>::reverse_iterator tti = tt.rbegin();
	for(int j=0; j<ntrial; j++)
	{		
		vector<IPoint>* pts= _buf->objAtTime(tti->second);
		tti++;

		int n = pts->size();
		if(n<m)
			continue;

		double e2 = 0.0;
		for(int z=0;z<n-1;z++)
			e2 += PDIST2(*(pts->at(z).p),*(pts->at(z+1).p));
		e2 /= n;
		//e2 = PDIST2(*(pts->at(0).p),*(pts->at(1).p));

		int cnt=0;
		while(true)
		{
			cnt++;
			int found=0; // number of objects found within the current eps

			// range query	(random point is the first one)
			int nloop = (n>nsample) ? nsample : n;
			for (int i=0; i<nloop; i++)
				if(PDIST2(*(pts->at(0).p),*(pts->at(i).p))<=e2)
					found++;

			if(found==0)
				e2 *= 1.9;
			else if (found>(int)m)
				e2  /= 2;
			else
			{
				double r = (double) n/_buf->num_obj();
				e2 *= r*r;
				break;
			}
			if (cnt>stop)
				break;
			//assert(cnt<stop);
		}
		es.insert(e2);
	}
	
	set<double>::iterator it = es.begin();	
	//int nloop = (int)es.size()/2; // median
	int nloop = (int)(es.size()*k/_buf->num_tstamp()); // T : T' = k : k'
	if(nloop>es.size()) nloop=es.size();
	for(int z=0; z<nloop; z++)
		it++;
	return sqrt(*it);
}

void ConvoyQuery::output( Convoy& c )
{
	if(_rst.empty())
	{
		_rst.push_back(c);
		return;
	}
	else
	{ 
		for(int j=0; j<_rst.size(); j++)
		{
			if( (c._te < _rst[j]._ts) || (_rst[j]._te < c._ts))
				continue;   

			set<int> is;
			set_intersection(_rst[j].begin(),_rst[j].end(), c.begin(),c.end(),
				insert_iterator<set<int>>(is,is.begin()));	
			if (is.size()==c.size())
			{
				_rst[j]._ts = MIN(_rst[j]._ts,c._ts);
				_rst[j]._te = MAX(_rst[j]._te,c._te);						
				return;
			}
		}
		_rst.push_back(c);
	}	
}
