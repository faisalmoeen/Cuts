
/*---------------------------------------------------------------------
DESCRIPTION: Coherent Moving Cluster (CMC)
perform a convoy query with the MC2 algorithm which
was introduced in "On discovering moving clusters
in spatio-temporal data"-SSTD05 paper. However, we
modify the algorithm to be able to find not moving
clusters but cv.
In convoy discovery, the theta value must be 1, thus
while clause in the original algorithm is not needed
here.
PARAMETERS
- int m : minimum number of objects
- int k : minimum lifetime
- double e : minimum closeness
AUTHOR: Hoyoung Jeung, 4 June 2007
TODO: do we need the continuousness time check? can it be pruned by
g->assigned ?
*---------------------------------------------------------------------*/


//#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "buf.h"
#include "func.h"
//#include "tmp.h"
#include "cmc.h"
#include "mc.h"
#include "cuts.h"

#include "sstream"
#include "string"


using namespace std;

int main()
{

	ifstream fo_Par("ConvoyP.txt");
	string Parameters[10];

	// For several experiments
	int ii = 0;
	for (std::string line; getline(fo_Par, line); )
	{
		// Current experiment with its set of parameters
		stringstream ActExp(line);
		//  Extract the parameters for each experiment
		int jj = 0;
		for (std::string par;getline(ActExp, par, ',');)
		{
			Parameters[jj] = par;
			cout << "Param: " << par << endl;
			jj++;
		}

		// -------------------- Once we have the parameters we set the variables

		Buffer buf;
		ofstream fo;
		double tol;
		int plen, tol_s, tol_e, tol_i, q;
		int		m = stoi(Parameters[0]); // ordered by truck, cattle, car, taxi
		int		k = stoi(Parameters[1]);
		double	e = stoi(Parameters[2]);
		int t = stoi(Parameters[3]);  // tolerance from manual
		int	l = stoi(Parameters[4]); // lambda from manual

									 // ----------------- Perform the experiment
		vector<char*> db;
		db.push_back("data/truck_min.txt");	 // N=267, T=10586 (max 528, avg 224), 59894 pts

		ii++;  // To count the number of experiments for the output file name
		fo.open("cmc" + to_string(ii) + ".txt", ios::out);
		Buffer b(db[0]); fo << db[0] << endl;
		CoherentMovingCluster cmc(&b); cmc.useIndex(false);
		cmc.discover(m, k, e);fo << "CMC "; cmc.summary(fo); fo << endl;
		CuTS c(&b); c.lambda(l);
		// Simple CuTS
		c.star(false);c.plus(false);c.simplify(t);c.discover(m, k, e); c.summary2(fo);fo << endl;
		// Simple CuTS +
		c.star(false);c.plus(true);c.simplify(t);c.discover(m, k, e); c.summary2(fo); fo << endl;
		// Simple CuTS *
		c.star(true);c.plus(false);c.simplify(t);c.discover(m, k, e); c.summary2(fo);fo << endl;
		fo << "-------------" << endl;
	}
	
	return 0;
}

/*#include <cstdio>


int main()
{
    printf("hello from ConvoyForLinux!\n");
    return 0;
}
*/