/*
#include "func.h"

func::func()
{
}
func::~func()
{
}
*/
//#include "stdafx.h"
#include "func.h"
#include <time.h>

Function::Function() {}
Function::~Function() {}

/*---------------------------------------------------------------------
DESCRIPTION: compute average velocity during from reftm
to (reftm-numRetro)
PARAMETERS
1. reftm: reference time
2. retro: retrospect-the number of time stamp in history from the
reference time
RETURN
1. double : velocity.
AUTHOR: Hoyoung Jeung, 2,nov, 2005
*---------------------------------------------------------------------*/
double	Function::velocity(Point *pts, int reftm, int retro)
{
	double dist = 0.0;
	for (int i = reftm; i>reftm - retro; i--)
	{
		dist += PDIST2(pts[i], pts[i - 1]);
	}
	return dist / retro;
}


/*---------------------------------------------------------------------
DESCRIPTION: drand48() functions return non-negative, double
precision, doubleing-point values, uniformly distributed over
the interval [0.0 , 1.0]. It originates from C language but
C++ doesn't have this function. That's why we define it.
FlockGenParams
1.void
RETURN
1.doubleing function number
AUTHOR: Hoyoung Jeung, 22, sept, 2005
NOTES:
*---------------------------------------------------------------------*/
double Function::drand48()
{
	return uniform(0.0, 1.0);
}


/************************************************************
** Generates a function number between _min and _max         **
** uniformly                                               **
By Yufei Tao
************************************************************/
double Function::uniform(double _min, double _max)
{
	int int_r = rand();
	int base = RAND_MAX - 1;
	double f_r = ((double)int_r) / base;
	return (_max - _min) * f_r + _min;
}

double Function::new_uniform(int _d_num)
{
	double base = 1;
	double sum = 0;
	for (int i = 0; i<_d_num; i++)
	{
		int digit = (int)uniform(0, 10);
		if (digit == 10) digit = 9;
		sum += base*digit;
		base *= 10;
	}
	return sum;
}

double Function::new_uniform(double _min, double _max)
{
	double ran_base = 9999999;
	double ran = new_uniform(7);
	return ran / ran_base*(_max - _min) + _min;
}



/*---------------------------------------------------------------------
DESCRIPTION: generate a random integer number which is not
duplicate to avoid
PARAMETERS
1. _min: minimum range
2. _max: maximum range
3. avoid: the result random number must not be this number
RETURN
1. int: a number in the range of [_min, _max] and not equals to avoid
AUTHOR: Hoyoung Jeung, 30 June 2006
*---------------------------------------------------------------------*/

int Function::uniform(int _min, int _max, int avoid)
{
	double rst;
	while (rst == avoid || rst<_min || rst>_max)
		rst = uniform((double)_min, (double)_max);

	return (int)rst;
}

int Function::uniform(int _min, int _max)
{
	return (int)uniform((double)_min, (double)_max);
}


/************************************************************
***  Given a mean and a standard deviation(sigma), gaussian       **
**   generates a normally distributed function number        **
**   Algorithm:  Polar Method, p.  104, Knuth, vol. 2      **
************************************************************/
double Function::gaussian(double mean, double sigma)
{
	double v1, v2;
	double s;
	double x;

	do
	{
		v1 = 2 * uniform(0, 1) - 1;
		v2 = 2 * uniform(0, 1) - 1;
		s = v1*v1 + v2*v2;
	} while (s >= 1.);

	x = v1 * sqrt(-2. * log(s) / s);

	/*  x is normally distributed with mean 0 and sigma 1.  */
	x = x * sigma + mean;

	return (x);
}


/*************************************************************/
/*  zipf generates a function number that follows Zipf         **
**  distribution and lies between x1 and x2.                 **
**  original code by Christos Faloutsos, 1995

**  The original node outputs discrete data only. The current**
**  function remedies the problem.			                 **
**  Modified by Yufei Tao (08/Dec/02)                         **
**************************************************************/
double Function::zipf(double x1, double x2, double p)
{
	double x;
	double i;
	double r, HsubV, sum;
	int V = 100;

	//double uniform();

	/* calculate the V-th harmonic number HsubV. WARNING: V>1 */
	HsubV = 0.0;
	for (i = 1; i <= V; i++)
		HsubV += 1.0 / pow((double)i, p);

	r = uniform(0., 1.)*HsubV;
	sum = 1.0; i = uniform(0, 1);
	while (sum<r) {
		//i++;  //commented by Yufei Tao
		i += uniform(1, 2);
		sum += 1.0 / pow((double)i, p);
	}

	/* i follows Zipf distribution and lies between 1 and V */

	/* x lies between 0. and 1. and then between x1 and x2 */
	x = ((double)i - 1.) / ((double)V - 1.);
	x = (x2 - x1) * x + x1;

	return(x);
}


/*---------------------------------------------------------------------
DESCRIPTION: count the number of lines
PARAMETERS
1. fname: a file name to count
RETURN
1. int: the number of lines
AUTHOR: Hoyoung Jeung, 19 Apr 2007
*---------------------------------------------------------------------*/
int Function::numLine(char *fname)
{
	int cnt = 0;
	char c;
	FILE *fp = fopen(fname, "r"); //C4996
	if (!fp)
	{
		printf("Function::numLine: cannot open the given file\n");
		exit(0);
	}

	do
	{
		c = fgetc(fp);
		if (c == '\n')
			cnt++;
	} while (c != EOF);

	fclose(fp);
	return cnt;
}
/*---------------------------------------------------------------------
DESCRIPTION: get the map extent from a given file
ConvoyGenParams
1. fn: input(source dataset) file name;
2. nline: # of the file line
3-6. containers of the max/min variables
RETURN
AUTHOR: Hoyoung Jeung, 27, Feb, 2006
NOTES: input MUST BE below form!
0 11.111 22.222
1 33.333 44.444
:
*---------------------------------------------------------------------*/
void Function::extent(char *fn, int nline,
	double *maxx, double *maxy, double *minx, double *miny)
{
	FILE *fip = fopen(fn, "r");
	if (!fip)
	{
		printf("Database.extent: could not read the file\n");
		return;
	}
	double x, y;
	int t;

	for (int i = 0; i<nline; i++)
	{
		fscanf(fip, "%d %lf %lf\n", &t, &x, &y);
		if (i == 0)
		{
			*maxx = *minx = x;
			*maxy = *miny = y;
		}
		*maxx = MAX(*maxx, x);
		*minx = MIN(*minx, x);
		*maxy = MAX(*maxy, y);
		*miny = MIN(*miny, y);
	}
}

void Function::progress(int now, int total)
{
	int ten = total / 10;
	if (now<ten - 1) return;
	else if (now == ten)   cout << "10% - done" << endl;
	else if (now == ten * 2) cout << "20% - done" << endl;
	else if (now == ten * 3) cout << "30% - done" << endl;
	else if (now == ten * 4) cout << "40% - done" << endl;
	else if (now == ten * 5) cout << "50% - done" << endl;
	else if (now == ten * 6) cout << "60% - done" << endl;
	else if (now == ten * 7) cout << "70% - done" << endl;
	else if (now == ten * 8) cout << "80% - done" << endl;
	else if (now == ten * 9) cout << "90% - done" << endl;
	else if (now == total) cout << "100% - done" << endl;

}

/*---------------------------------------------------------------------
DESCRIPTION:
- Trajectory& trj :
AUTHOR: Hoyoung Jeung, 5/11/2007
NOTE:
*---------------------------------------------------------------------*/
void Function::sort(vector<Point>* trj)
{
	Point temp;
	for (int i = 0; i < trj->size() - 1; i++)
		for (int j = i + 1; j < trj->size(); j++)
			if (trj->at(i).t>trj->at(j).t)
			{
				temp = trj->at(i);    //swapping entire struct
				trj->at(i) = trj->at(j);
				trj->at(j) = temp;
			}
}

