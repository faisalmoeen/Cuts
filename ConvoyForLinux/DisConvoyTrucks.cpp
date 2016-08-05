
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


#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "buf.h"
#include "func.h"
#include "tmp.h"
#include "cmc.h"
#include "mc.h"
#include "cuts.h"

#include "sstream"
#include "string"



using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{

	ifstream fo_Par("ConvoyP.txt");
	string Parameters[10];

	// For several experiments
	int ii=0;
	for( std::string line; getline( fo_Par, line ); )
		{
			// Current experiment with its set of parameters
			stringstream ActExp(line);
			//  Extract the parameters for each experiment
			int jj=0;
			for (std::string par;getline(ActExp, par,',');)
			{
				Parameters[jj]=par;
				cout <<"Param: " << par<<endl;
				jj++;
			}
	
			// -------------------- Once we have the parameters we set the variables

				Buffer buf;
				ofstream fo;
				double tol;
				int plen,tol_s,tol_e,tol_i, q;
				int		m = stoi(Parameters[0]); // ordered by truck, cattle, car, taxi
				int		k = stoi(Parameters[1]);
				double	e = stoi(Parameters[2]);
				int t = stoi(Parameters[3]);  // tolerance from manual
				int	l = stoi(Parameters[4]); // lambda from manual

				// ----------------- Perform the experiment
				vector<char*> db;
				db.push_back("data/truck_min.txt");	 // N=267, T=10586 (max 528, avg 224), 59894 pts

				ii++;  // To count the number of experiments for the output file name
				fo.open("cmc" + to_string(ii) + ".txt",ios::out);
				Buffer b(db[0]); fo << db[0] << endl;	
				CoherentMovingCluster cmc(&b); cmc.useIndex(false);
				cmc.discover(m,k,e);fo << "CMC "; cmc.summary(fo); fo << endl;
				CuTS c(&b); c.lambda(l);
				// Simple CuTS
				c.star(false);c.plus(false);c.simplify(t);c.discover(m,k,e); c.summary2(fo);fo << endl;
				// Simple CuTS +
				c.star(false);c.plus(true);c.simplify(t);c.discover(m,k,e); c.summary2(fo); fo << endl;
				// Simple CuTS *
				c.star(true);c.plus(false);c.simplify(t);c.discover(m,k,e); c.summary2(fo);fo << endl;
				fo << "-------------" << endl;
	}


	/*
	Buffer buf;
	ofstream fo;
	double tol;
	int plen,tol_s,tol_e,tol_i, q;

	//	The first position of this vectors has the information for trucks.

	int		m[] = {6  ,2  ,3  ,3}; // ordered by truck, cattle, car, taxi
	int		k[] = {180,180,180,180};
	double	e[] = {200 ,300,80 ,40};

//	double  t[] = {5.9 ,274.2,63.4 ,31.5};  // tolerance : from automation
//	int	l[] = {4  ,36  ,24 ,4}; // lambda from automation

//	int t[] = {8  ,84 ,7  ,26};  // tolerance from manual
	//int	l[] = {8  ,30 ,30 ,10}; // lambda from manual

	// Other parameters
	int t[] = {8  ,84 ,7  ,26};  // tolerance from manual
	int	l[] = {8  ,30 ,30 ,10}; // lambda from manual


	//fo.open("truck.txt",ios::out); // Out for output operations

	cout << "CMC Clustering";
	vector<char*> db;
	db.push_back("data/truck_min.txt");	 // N=267, T=10586 (max 528, avg 224), 59894 pts

	cout<<"----- ojo - ";
	cout<< db[0] << " / fin ";
	//cout<< db[1] << " / fin2 ";

	int z=15;

	if (z==15) //----- exp : cmc vs cuts  -----
	{
		fo.open("cmc.txt",ios::out);
		for(int q=0; q<=0; q++){
			Buffer b(db[q]); fo << db[q] << endl;	
			CoherentMovingCluster cmc(&b); cmc.useIndex(false);
			cmc.discover(m[q],k[q],e[q]);fo << "CMC "; cmc.summary(fo); fo << endl;
			CuTS c(&b); c.lambda(l[q]);
			c.star(false);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl;
			c.star(false);c.plus(true);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo); fo << endl;
			c.star(true);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl;
			fo << "-------------" << endl;
		}
	}
	*/

	/*
  if (z==4) //----- test : cmc test -----
	{	
		cout << "_____________________Llegamos Will: ";
		q = 0;
		buf.open(db[q]); 
		CoherentMovingCluster cmc(&buf);
		//cmc.useIndex(false);cmc.discover(m[q],k[q],e[q]);cmc.write("data/test/cmc_noidx.txt");	
		cmc.useIndex(true);cmc.discover(m[q],k[q],e[q]);cmc.write("data/test/cmc.txt");	
		vector<Convoy>* r = cmc.getResults();
		int sum1=0;
		int sum2=0;
		for(int i=0; i<r->size(); i++)
		{
			sum1 +=r->at(i).size();
			sum2 +=r->at(i).lifetime();
		}
		cout << "avg convoy objects=" << (double)sum1/r->size() 
			 << " , avg lifetime="<< (double)sum2/r->size() << endl;
	}			



if (z==-1) //----- test : cmc test -----
	{
		fo.open("truck.txt",ios::out);
		for(int q=0; q<=0; q++){	
			Buffer b(db[q]);fo << db[q] << endl; cout << db[q] << endl;
			CoherentMovingCluster cmc(&b);
			cmc.useIndex(true);
			cmc.discover(m[q],k[q],e[q]);//cmc.write("data/test/cmc.txt");			
			CuTS c(&b); c.lambda(b.num_tstamp());	
			c.star(false);c.plus(false);c.simplify(0);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 			
			c.star(false);c.plus(true);c.simplify(0);c.discover(m[q],k[q],e[q]); c.summary2(fo); fo << endl;  
			c.star(true);c.plus(false);c.simplify(0);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 
			fo << "CMC : "; cmc.summary(fo);			
		}
	}	
	*/
	


/*	if (z==-1) //----- test : cmc test -----
	{
		fo.open("truck.txt",ios::out);
		for(int q=0; q<=0; q++){	
			Buffer b(db[q]);fo << db[q] << endl; cout << db[q] << endl;
//			CoherentMovingCluster cmc(&b);
//			cmc.useIndex(true);
//			cmc.discover(m[q],k[q],e[q]);//cmc.write("data/test/cmc.txt");			
			cout << "entra CUTS: "<< endl;
			CuTS c(&b); 
			c.lambda(b.num_tstamp());	
			c.star(false);
			c.plus(false);
			c.simplify(0);
					cout << "Fin Simplify " << endl;
			c.discover(m[q],k[q],e[q]); 
					cout << "Fin Discover " << endl;
			c.summary2(fo);fo << endl; 			
			cout << "Fin CUTS 1: "<< endl;

			c.star(false);
			c.plus(true);
			c.simplify(0);
			c.discover(m[q],k[q],e[q]);
			c.summary2(fo); 
			fo << endl;  
			
			cout << "Fin CUTS 2: "<< endl;
			c.star(true);
			c.plus(false);
			c.simplify(0);
			c.discover(m[q],k[q],e[q]);
			c.summary2(fo);
			fo << endl; 
//			fo << "CMC : "; cmc.summary(fo);			
		}
	}	
	*/
	/*
	// OK - Execute simplification with DP - DP+ and DP*
	if (z==13) //----- exp : simplification -----
	{
		fo.open("dp.txt",ios::out); 
		for(int q=0; q<=0; q++){
			Buffer b(db[q]); fo << db[q] << endl; cout << db[q] << endl;
			CuTS c(&b);	c.useIndex(false);  c.lambda(l[q]);
			int t,r;
			for(int s=10; s<=50; s+=10){fo << "tol=" << s << endl;
				c.star(false);c.plus(false); t=c.simplify(s); r=c.reduction();
				fo<<"DP  : "<< t <<" ms, "<< r << " %\n";
				c.star(false);c.plus(true);t=c.simplify(s); r=c.reduction();
				fo<<"DP+  : "<< t <<" ms, "<< r << " %\n";
				c.star(true);c.plus(false);t=c.simplify(s); r=c.reduction();
				fo<<"DP*  : "<< t <<" ms, "<< r << " %\n";
			}
			fo << "-------------" << endl;
		}
	}
	
	
	if (z==14) //----- exp : mc2 vs cmc -----
	{
		fo.open("mc2.txt",ios::out);
		for(int q=0; q<=0; q++){	
			Buffer b(db[q]); fo << endl << db[q] << endl;	
			CoherentMovingCluster cmc(&b); cmc.useIndex(true);
			cmc.discover(m[q],k[q],e[q]); fo << "CMC : "; cmc.summary(fo);
			vector<Convoy> cmcrst = *(cmc.getResults());		
			for(double theta=0.4; theta<=1.0; theta+=0.2){				 				
				cout << "MC2 : theta=" << theta << endl; 				
				cmc.mc2(theta,m[q],k[q]); cmc.validate(cmcrst);
				fo << "MC2 (" << theta << ") : "; cmc.summary(fo);
			}
		}
	}
	*/
/*
	if (z==15) //----- exp : cmc vs cuts  -----
	{
		fo.open("cmc.txt",ios::out);
		for(int q=0; q<=3; q++){	
			Buffer b(db[q]); fo << db[q] << endl;	
			CoherentMovingCluster cmc(&b); cmc.useIndex(false);
			cmc.discover(m[q],k[q],e[q]);fo << "CMC "; cmc.summary(fo); fo << endl;
			CuTS c(&b); c.lambda(l[q]);	
			c.star(false);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 			
			c.star(false);c.plus(true);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo); fo << endl;  
			c.star(true);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 
			fo << "-------------" << endl;
		}
	}
	if (z==16) //----- exp : actual tolerance -----
	{
		fo.open("act_tol.txt",ios::out); 
		for(int q=0; q<=3; q++){	
			Buffer b(db[q]); fo << db[q] << endl;
			CuTS c(&b);	c.useIndex(false);  c.lambda(l[q]);
				for(int at=0; at<=1; at++){
				c.actualTol(at);
				c.star(false);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 			
				c.star(false);c.plus(true);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo); fo << endl;  
				c.star(true);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 
				fo << endl;		
			}
			fo << "-------------" << endl;
		}
	}
	if (z==17) //----- exp : lambda -----
	{
		fo.open("lambda.txt",ios::out); 
		for(int q=0; q<=3; q++){	
			Buffer b(db[q]); fo << db[q] << endl;
			CuTS c(&b);	c.useIndex(false); 
			if(q==0 || q==3)
				for(int s=5; s<=20; s+=5){
					c.lambda(s);
					c.star(false);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 			
					c.star(false);c.plus(true);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo); fo << endl;  
					c.star(true);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 
					fo << endl;	}
			else 
				for(int s=10; s<=70; s+=20){
					c.lambda(s);
					c.star(false);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 			
					c.star(false);c.plus(true);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo); fo << endl;  
					c.star(true);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 
					fo << endl; }
			fo << "-------------" << endl;
		}
	}

	if (z==18) //----- exp : tol -----
	{
		fo.open("tol.txt",ios::out); 
		for(int q=2; q<=3; q++){	
			Buffer b(db[q]); fo << db[q] << endl;
			CuTS c(&b);	c.useIndex(false);  c.lambda(l[q]);
			if(q==0)
				for(int s=5; s<=20; s+=5){
					c.star(false);c.plus(false);c.simplify(s);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 			
					c.star(false);c.plus(true);c.simplify(s);c.discover(m[q],k[q],e[q]); c.summary2(fo); fo << endl;  
					c.star(true);c.plus(false);c.simplify(s);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 
					fo << endl; }
			else 
				for(int s=10; s<=250; s+=70){
					c.star(false);c.plus(false);c.simplify(s);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 			
					c.star(false);c.plus(true);c.simplify(s);c.discover(m[q],k[q],e[q]); c.summary2(fo); fo << endl;  
					c.star(true);c.plus(false);c.simplify(s);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 
					fo << endl; }
			fo << "-------------" << endl;
		}
	}
	if (z==74) //----- test : e -----
	{
		fo.open("e.txt",ios::out); 
		for(int q=0; q<=3; q++){	
			Buffer b(db[q]); fo << db[q] << endl;
			//CoherentMovingCluster cmc(&b); cmc.useIndex(false);
			//cmc.discover(m[q],k[q],e[q]);fo << "CMC "; cmc.summary(fo); fo << endl;
			CuTS c(&b);	c.useIndex(false);  c.lambda(l[q]);
			for(double s=1.0; s<=1.9; s+=0.3){
				c.star(false);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]*s); c.summary2(fo);fo << endl; 			
				c.star(false);c.plus(true);c.simplify(t[q]);c.discover(m[q],k[q],e[q]*s); c.summary2(fo); fo << endl;  
				c.star(true);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]*s); c.summary2(fo);fo << endl; 
				fo << endl; 
			}
			fo << "-------------" << endl;
		}
	}
	if (z==2) //----- test : cuts -----
	{	
		fo.open("car_week.txt",ios::out);
		q = 2;	buf.open(db[q]);  fo << db[q] << endl;			 
		CoherentMovingCluster cmc(&buf); 
		cmc.discover(m[q],k[q],e[q]);fo << "CMC "; cmc.summary(fo); fo << endl;
		cmc.useIndex(true);cmc.discover(m[q],k[q],e[q]);fo << "CMC "; cmc.summary(fo); fo << endl;
		CuTS c(&buf); c.lambda(l[q]);
		c.star(false);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 			
		c.star(false);c.plus(true);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo); fo << endl;  
		c.star(true);c.plus(false);c.simplify(t[q]);c.discover(m[q],k[q],e[q]); c.summary2(fo);fo << endl; 
	}
	if (z==90) //----- test : batch  -----
	{
		fo.open("taxi1000_batch.txt",ios::out);
		for(int q=3; q<=3; q++){	
			Buffer b(db[q]); fo << db[q] << endl;	
			CuTS c(&b);	c.useIndex(false);  c.lambda(l[q]);
			for(int s=0; s<=2; s++){
				c.batch(s);	c.simplify(t[q]);
				c.discover(m[q],k[q],e[q]); c.summary2(fo);	fo << endl; 
			}
			fo << "-------------" << endl;
		}
	}
	if (z==-3) //----- test : delta and lambda computation -----
	{	
		int		N[] = {267  ,13     ,183  ,8757}; // ordered by truck, cattle, car, taxi
		int		T[] = {10586,175636 ,8757 ,965};
		double	A[] = {224  ,175636 ,451  ,82};
		double  D[] = {59894 ,2283268,82590 ,41144}; 
		int		K[] = {3 ,180,3 ,3}; 
		int reduc = 50;
		for(int q=0; q<=3; q++){	
			//cout << (int) (2 * T[q] * N[q]/D[q] * 100/reduc)
			int xx = (1-k[q]/A[q])*10;
			xx = (xx<=0) ? 1 : xx;
			cout << (int) (2 * xx * 100/reduc)
				<< endl;
		}
	}
	if (z==24) //----- test : delta and lambda computation -----
	{	
		for(int q=0; q<=3; q++){	
			Buffer b(db[q]);fo << db[q] << endl << endl; cout << db[q] << endl << endl;
			CuTS c(&b);
			double tol = c.compute_tol(e[q]);
			c.simplify(t[q]);
			double lam = c.computeLambda(m[q],k[q],c.reduction());
			//double lam = c.computeLambda(m[q],k[q],50);
			cout << "tol=" << tol << ", lambda=" << lam << endl;
			//cout << ", lambda=" << lam << endl;
			c.summary2(fo);fo << endl; 
		}
	}	
	if (z==26) //----- test : index -----
	{
		fo.open("taxi_idx.txt",ios::out);
		for(int q=3; q<=3; q++){	
			Buffer b(db[q]); fo << db[q] << endl;	
			CuTS c(&b);
			c.lambda(8);
			for(int s=0; s<=2; s++){
				c.useIndex(s); c.lambda(l[q]); c.simplify(t[q]); 
				c.discover(m[q],k[q],e[q]); c.summary2(fo);	
				fo << endl; 
			}
			fo << "-------------" << endl;
		}
	}

	if (z==4) //----- test : cmc test -----
	{	
		cout << "_____________________Llegamos Will: ";
		q = 2;
		buf.open(db[q]); 
		CoherentMovingCluster cmc(&buf);
		//cmc.useIndex(false);cmc.discover(m[q],k[q],e[q]);cmc.write("data/test/cmc_noidx.txt");	
		cmc.useIndex(true);cmc.discover(m[q],k[q],e[q]);cmc.write("data/test/cmc.txt");	
		vector<Convoy>* r = cmc.getResults();
		int sum1=0;
		int sum2=0;
		for(int i=0; i<r->size(); i++)
		{
			sum1 +=r->at(i).size();
			sum2 +=r->at(i).lifetime();
		}
		cout << "avg convoy objects=" << (double)sum1/r->size() 
			 << " , avg lifetime="<< (double)sum2/r->size() << endl;
	}			

	if (z==3) // ------ test : distance functions
	{
		PolyLine l1(1),l2(2); 
		l1.push_back(&Point(0,1,1));l1.push_back(&Point(1,4,3));l2.push_back(&Point(0,2,3));l2.push_back(&Point(1,3,0));
		//l1.push_back(&Point(0,2,2));l1.push_back(&Point(1,6,2));l2.push_back(&Point(0,2,1));l2.push_back(&Point(1,6,2));
		//l1.push_back(&Point(0,4,4));l1.push_back(&Point(1,6,6));l2.push_back(&Point(0,4,6));l2.push_back(&Point(1,5,5));
		//l1.push_back(&Point(0,1,1));l1.push_back(&Point(1,3,1));l2.push_back(&Point(0,1,2));l2.push_back(&Point(1,3,2));
		//l1.push_back(&Point(0,5,1));l1.push_back(&Point(1,5,3));l2.push_back(&Point(0,6,1));l2.push_back(&Point(1,6,3));
		//l1.push_back(&Point(0,2,6));l1.push_back(&Point(1,4,8));l2.push_back(&Point(0,2,8));l2.push_back(&Point(1,4,6));
		//l1.push_back(&Point(0,1,1));l1.push_back(&Point(1,4,3));l1.push_back(&Point(2,5,1));l1.push_back(&Point(3,1,0));
		//l2.push_back(&Point(0,2,3));l2.push_back(&Point(1,3,0));l2.push_back(&Point(2,2,1));l2.push_back(&Point(3,2,0));
		double d,d_st,d_box;
		l1.distLL(l2,&d); l1.distLL_ST(l2,&d_st);l1.ext.min_dist(l2.ext,&d_box);
		cout << "d=" << d << ", d_st=" << d_st << ", d_box=" << d_box << endl;
	}
	if (z==11) // ------ test : mbr intersection
	{
		//Point p1(0,2,5), p2(0,1,5);	MBR m(1,4,3,5);cout << "intersect = " << m.intersect(p1,p2) <<endl;
		//PolyLine l(1);l.push_back(&Point(0,1,1));l.push_back(&Point(1,1,3));l.push_back(&Point(1,2,4));
		//l.push_back(&Point(1,2,4.5));l.push_back(&Point(1,2,5.5));l.push_back(&Point(1,3,6));
		//MBR m(1,4,3,5);
		//for(int z=0; z<l.size()-1;z++)
		//	cout << "intersect = " << m.intersect(*l[z],*l[z+1]) <<endl;
	}

	if (z==100) //----- data processing -----
	{	
		//Buffer::truck("data/raw/trucks_raw.txt","data/truck_min.txt");		
		//Buffer::timestat("data/truck_min.txt","data/truck_min_stat.txt");
		//Buffer::truck("data/raw/trucks_raw.txt","data/truck_30sec.txt");		
		//Buffer::timestat("data/truck_30sec.txt","data/truck_30sec_stat.txt");
		
		//Buffer::beijing("data/raw/beijing14_raw.txt","data/taxi100_min.txt",100);
		//Buffer::timestat("data/taxi100_min.txt","data/taxi100_min_stat.txt");
		Buffer::beijing("data/raw/beijing14_raw.txt","data/taxi500_min.txt",500);
		Buffer::timestat("data/taxi500_min.txt","data/taxi500_min_stat.txt");
		//Buffer::beijing("data/raw/beijing14_raw.txt","data/taxi1000.txt",1000);
		//Buffer::timestat("data/taxi1000.txt","data/taxi1000_stat.txt");
		
		//Buffer::akta2("data/raw/akta","data/car_day_min.txt",2,1,15,15); // one day
		//Buffer::timestat("data/car_day_min.txt","data/car_day_min_stat.txt");	
		//Buffer::akta2("data/raw/akta","data/car_week_min.txt",1,12,17,23); // one week - min
		//Buffer::timestat("data/car_week_min.txt","data/car_week_min_stat.txt");
		//Buffer::akta2("data/raw/akta","data/car_sec_week.txt",2,1,25,31); // one week - sec
		//Buffer::timestat("data/car_sec_week.txt","data/car_week_sec_stat.txt");
		//Buffer::akta2("data/raw/akta","data/car_month_min.txt",1,10,1,30); // one month - min
		//Buffer::timestat("data/car_month_min.txt","data/car_month_min_stat.txt");	
		//Buffer::extractSameTime("data/truck.txt","data/truck_.txt");	
	}

fo.close();
	*/


	return 0;
}
