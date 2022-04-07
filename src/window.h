#ifndef WINDOW_H_
#define WINDOW_H_
#include "base.h"
#include "util.h"
#include "kvp.h"
#include "vec.h"
#include "mat.h"
#include "polytope.h"
#include "hitandrun.h"
#include "optim.h"
#include "redund.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#define RADIUS .75,
using namespace std;
template <typename T>
class window
{
public:
	//constructor
	mat<T> lb;        // vector of size_x lower bounds contained in a vector over T frames, contained in a vector of N trajectories
	mat<T> ub;        // vector of size_x upper bounds contained in a vector over T frames, contained in a vector of N trajectories
	mat<T> traj;      // vector of size_x states contained in a vector over T frames, contained in a vector of N trajectories
	vec<T> START;
	mat<T> skew;
	vec<T> current;
	std::string path;
	std::ostringstream ss;
	std::size_t F;//number of time frames
	T DT = 0;
	T TCD = 0;
	T TCA = 0;
	std::size_t size_x; //dim of A, b
	std::size_t size_y;
	std::size_t HARN;
	vec<bool> insides;
	window(std::size_t Fin,T DTin,T TCAin,T TCDin, const std::string& pathin)//,std::size_t Nin
	{
		F = Fin;
		DT = DTin;
		TCA= TCAin;
		TCD = TCDin;
		path = pathin;
		mat<T> H(path+"000.H");
		size_y = H.size_y;
		size_x = H.size_x-1;
		HARN = std::size_t(pow(size_y,2));
		//cout << "size:"<<size_y << "," << size_x << endl;
		lb   = mat<T>(F,size_x);
		ub   = mat<T>(F,size_x);
		traj = mat<T>(F,size_x);
		skew = mat<T>(path+"deltas.csv");
		//cout << "skew:" << endl;
    //skew.print();
		insides = vec<bool>(F);
		START = vec<T>(size_x);
		current = vec<T>(size_x);
	};
	window()
	{
		skew=0;
		F = 1;
		DT = 1;
		TCA = 1;
		TCD = 1;
		path = "";
		START = vec<T>(1);
		current = vec<T>(1);
		skew = mat<T>(1,1);
		size_y = 1;
		size_x = 1;
		HARN = std::size_t(pow(size_y,2));
		lb   = mat<T>(1,1);
		ub   = mat<T>(1,1);
		traj = mat<T>(1,1);
	};
	window(const window<T>& WND)
	{
		skew=WND.skew;
		F = WND.F;
		DT = WND.DT;
		TCA = WND.TCA;
		TCD = WND.TCD;
		path = WND.path;
		size_y = WND.size_y;
		size_x = WND.size_x;
		lb   = WND.lb;
		ub   = WND.ub;
		traj = WND.traj;
		current = WND.current;
	};
	window<T>& operator=(const window<T>& WND)
	{
		skew=WND.skew;
		F = WND.F;
		DT = WND.DT;
		TCA = WND.TCA;
		TCD = WND.TCD;
		START = WND.START;
		path = WND.path;
		size_y = WND.size_y;
		size_x = WND.size_x;
		HARN = std::size_t(pow(size_y,2));
		lb   = WND.lb;
		ub   = WND.ub;
		traj = WND.traj;
		current = WND.current;
		return (*this);
	};
	//destructor
	~window()
	{
	};

	void optim_one(polytope<T>& P,std::size_t x,T v)
	{
		hitandrun<double> HAR(P, current, 1, 1);
		vec<T> obj(P.H.size_x);
		obj.set_all(0); obj(x+1) = v;
		vec<T> opt = P.H.lpsolve(obj);
		for (std::size_t X = 0; X < P.A.size_x; X++){
			current(X) = opt(X+1);
		}
	};

	void optim_all(polytope<T>& P,T v)
	{
		hitandrun<double> HAR(P, current, 1, 1);
		vec<T> obj(P.H.size_x);
		obj.set_all(v); obj(0) = 0;
		vec<T> opt = P.H.lpsolve(obj);
		for (std::size_t X = 0; X < P.A.size_x; X++){
			current(X) = opt(X+1);
		}
	};

	void walk(polytope<T>& P,std::size_t tf)
	{
	  hitandrun<double> HAR(P, current, 1, 1);
	  HAR.analyticalCenter();
	  //HAR.analyticalCenterNewtons();
	  //HAR.chebyshevCenter();
	  vec<T> bsk(P.A.size_x);
	  T pi = 2*acos(0.0);
	  T z =.67; //1.28;
	  current = HAR.current;
	  int N = (P.A.size_x) * (P.A.size_y);
		int n_steps= N;
	  //let z be 1.5 for easy hand math, 1.7 for 3 sds
	  n_steps = N;
	  if (tf>0){
		n_steps = N;
			for (std::size_t X = 0; X < P.A.size_x; X++){
				T dat = (1 - traj(tf-1,X))/TCA;           //excitation u 1
				T ddt = (0 - traj(tf-1,X))/TCD;           // excitation u 0
				T Q1 = traj(tf-1,X) + z*DT*ddt/2; // "estimate" the iqr1 by just median the lower range
				T Q2 = traj(tf-1,X);
				T Q3 = traj(tf-1,X) + z*DT*dat/2; // "estimate" the iqr3 by just median the upper range
				T bse = (Q3 + Q1 - 2*Q2)/(Q3-Q1); // skewness "estimate" using just range and mean, absolutely devilish!
				// T lb = traj(tf-1,X) + DT*ddt;
				// T ub = traj(tf-1,X) + DT*dat;
				// T sd = (ub - lb);
				//T bse = 3*(current(X) - traj(tf-1,X))/sd;
				if (bse > .99){
					bse = .99;
				}else if (bse < -.99){
					bse = -.99;
				}
				T g23 = pow(abs(bse),.6667);
				T dalph = sqrt((pi/2)*(g23/(g23 + pow(((4 - pi)/2),.6667)))); // MLE estimate of the shape factor
				//cout << "bse:" << bse << endl;
				//cout << "dalph:" << dalph << endl;
				T den = sqrt(1 - pow(dalph,2));
				//cout << "den:" << den << endl;
				bsk(X) = dalph/ (den);
				//cout << "bsk:" << bsk(X) << endl;
				if (bse > 0){ //fix the sign of skew (dalph is strictly positive)
					bsk(X) = -bsk(X);
				}
				bsk(X) =bsk(X) + skew(tf,X);
			//	cout << "bsk:" << bsk(X) << endl;
			}
		}



	  int DIM = int(P.A.size_x);
	  int NCON = int(P.b.size);
	  Eigen::VectorXd initialization(DIM);
	  Eigen::VectorXd sk(DIM);
	  Eigen::VectorXd median(DIM);
	  Eigen::MatrixXd A(NCON,DIM);
	  Eigen::VectorXd b(NCON);
		sk.setZero();
	  initialization.setZero();
	  median.setZero();
	  A.setZero();
	  b.setZero();
	  for (std::size_t Y = 0; Y < P.A.size_y; Y++){
	    for (std::size_t X = 0; X < P.A.size_x; X++){
	      A(Y,X) = P.A(Y,X);
	    }
	    b(Y) = P.b(Y);
	  }
	  if (tf >0){
		for (std::size_t X = 0; X < P.A.size_x; X++){
		    initialization(X) = current(X);
				sk(X) = bsk(X);
				median(X) = traj(tf-1,X);
		  }
		}else{
			for (std::size_t X = 0; X < P.A.size_x; X++){
			    initialization(X) = current(X);
					sk(X) = 0;
					median(X) = HAR.current(X);
			  }
		}

	  T radius =1;
	  Eigen::MatrixXd vo = generateDikinWalkSamples(initialization, A, b, radius, n_steps,sk, median);//
	  Eigen::VectorXd curr = vo.col(vo.cols()-1);
	  for (std::size_t X = 0; X < P.A.size_x; X++){
	    current(X) = curr(X);
	  }
		//writeToCSVfile("test.csv",vo.transpose());
	};

	mat<T> make_IDSYS(std::size_t frame)
	{
		mat<T> U(size_x,size_x+1);
		mat<T> L(size_x,size_x+1);
		ss.str("");
		ss << std::setfill('0') << std::setw(3) << frame;
		mat<T> IDSYS = (path+ss.str()+".H");
		ss.str("");
		for (std::size_t X = 0; X < L.size_y; X++){
			L(X,0,-lb(frame-1,X));
			L(X,X+1,1);
			U(X,0,ub(frame-1,X));
			U(X,X+1,-1);
		}
		std::size_t constrain_y = IDSYS.size_y;
		IDSYS = (IDSYS.concat(L).concat(U));
		mat<T> temp = IDSYS.redund();
		//temp.print();
		//fudge until the problem is solve-able
		int c;
		T this_fudge = FUDGE;
		while ((temp.size_x == 1) && (temp.size_y == 1)){
		for (std::size_t Y = 0; Y < constrain_y; Y++){//2*L.size_y
			IDSYS(Y,0)+=this_fudge;
		}
		temp = IDSYS.redund();
		this_fudge*=1.5;
		c++;
		if (c>MAXITER){
			cerr << "failed on redund." <<endl;
			exit(1);
		}
		}
		IDSYS = temp;
		return IDSYS;
	};

	mat<T> start_IDSYS()
	{
		mat<T> U(size_x,size_x+1);
		mat<T> L(size_x,size_x+1);
		mat<T> IDSYS = (path+"000.H");
		ss.str("");
		for (std::size_t X = 0; X < L.size_y; X++){
			L(X,0,-0);
			L(X,X+1,1);
			U(X,0,1);
			U(X,X+1,-1);
		}
		std::size_t constrain_y = IDSYS.size_y;
		IDSYS = (IDSYS.concat(L).concat(U));
		mat<T> temp = IDSYS.redund();
		//temp.print();
		//fudge until the problem is solve-able
		int c;
		T this_fudge = FUDGE;
		while ((temp.size_x == 1) && (temp.size_y == 1)){
			for (std::size_t Y = 0; Y < constrain_y; Y++){//2*L.size_y
				IDSYS(Y,0)+=this_fudge;
			}
		temp = IDSYS.redund();
		this_fudge*=1.5;
		c++;
		if (c>MAXITER){
			cerr << "failed on redund." <<endl;
			exit(1);
		}
		}
		IDSYS = temp;
		return IDSYS;
	};

	//use a good supplied guess to seed random (CMC does our job!). on row 1 of traj and compute ub lb frame 0
	void seed_trajectory()
	{
		mat<T> IDSYS = start_IDSYS();
		polytope<T> P(IDSYS);
		walk(P,0);
		insides(0) =  P.check_inside(current);
		//////////////////////////////////////////////////////////////////////////////
		// a = self.tau_activation*(.5 + 1.5*this_activation_set)
		// d = self.tau_deactivation/(.5 + 1.5*this_activation_set)
		// activation_lb = this_activation_set+delta_time*np.divide(0-this_activation_set,d)
		// activation_lb[activation_lb<0] = 0
		// activation_ub = this_activation_set+delta_time*np.divide(1-this_activation_set,a)
		// activation_ub[activation_ub>1] = 1
		//////////////////////////////////////////////////////////////////////////////
		for (std::size_t X = 0; X < size_x; X++){
			traj(0,X) = current(X);
			T dat = (1 - traj(0,X))/TCA;           //excitation u 1
			T ddt = (0 - traj(0,X))/TCD;           // excitation u 0
			ub(0,X) = traj(0,X) + DT*dat;
			lb(0,X) = traj(0,X) + DT*ddt;
			if (ub(0,X) > 1) ub(0,X) = 1; //set ubs greater than 1 to 1
			if (lb(0,X) < 0) lb(0,X) = 0; //set lbs less than 0 to 0
		}
	};

	void iter(std::size_t frame)
	{
		mat<T> IDSYS = make_IDSYS(frame); //use the previous frame
		polytope<T> P(IDSYS);
		walk(P,frame);
		insides(frame) =  P.check_inside(current);
		for (std::size_t X = 0; X < size_x; X++){
			traj(frame,X) = current(X);
			T dat = (1 - traj(frame,X))/TCA;           //excitation u 1
			T ddt = (0 - traj(frame,X))/TCD;           // excitation u 0
			ub(frame,X) = traj(frame,X) + DT*dat;
			lb(frame,X) = traj(frame,X) + DT*ddt;
			if (ub(frame,X) > 1) ub(frame,X) = 1; //set ubs greater than 1 to 1
			if (lb(frame,X) < 0) lb(frame,X) = 0; //set lbs less than 0 to 0
		}
	};

	//use a good supplied guess to seed random (CMC does our job!). on row 1 of traj and compute ub lb frame 0
	void seed_trajectory_optim_one(std::size_t m,T v)
	{
		mat<T> IDSYS = start_IDSYS();
		polytope<T> P(IDSYS);
		optim_one(P,m,v);
		for (std::size_t X = 0; X < size_x; X++){
			traj(0,X) = current(X);
			T dat = (1 - traj(0,X))/TCA;           //excitation u 1
			T ddt = (0 - traj(0,X))/TCD;           // excitation u 0
			ub(0,X) = traj(0,X) + DT*dat;
			lb(0,X) = traj(0,X) + DT*ddt;
			if (ub(0,X) > 1) ub(0,X) = 1; //set ubs greater than 1 to 1
			if (lb(0,X) < 0) lb(0,X) = 0; //set lbs less than 0 to 0
		}
	};

	void iter_optim_one(std::size_t frame,std::size_t m,T v)
	{
		mat<T> IDSYS = make_IDSYS(frame); //use the previous frame
		polytope<T> P(IDSYS);
		optim_one(P,m,v);
		for (std::size_t X = 0; X < size_x; X++){
			traj(frame,X) = current(X);
			T dat = (1 - traj(frame,X))/TCA;           //excitation u 1
			T ddt = (0 - traj(frame,X))/TCD;           // excitation u 0
			ub(frame,X) = traj(frame,X) + DT*dat;
			lb(frame,X) = traj(frame,X) + DT*ddt;
			if (ub(frame,X) > 1) ub(frame,X) = 1; //set ubs greater than 1 to 1
			if (lb(frame,X) < 0) lb(frame,X) = 0; //set lbs less than 0 to 0
		}
	};

	void seed_trajectory_optim_all(T v)
	{
		mat<T> IDSYS = start_IDSYS();
		polytope<T> P(IDSYS);
		optim_all(P,v);
		for (std::size_t X = 0; X < size_x; X++){
			traj(0,X) = current(X);
			T dat = (1 - traj(0,X))/TCA;           //excitation u 1
			T ddt = (0 - traj(0,X))/TCD;           // excitation u 0
			ub(0,X) = traj(0,X) + DT*dat;
			lb(0,X) = traj(0,X) + DT*ddt;
			if (ub(0,X) > 1) ub(0,X) = 1; //set ubs greater than 1 to 1
			if (lb(0,X) < 0) lb(0,X) = 0; //set lbs less than 0 to 0
		}
	};

	void iter_optim_all(std::size_t frame,T v)
	{
		mat<T> IDSYS = make_IDSYS(frame); //use the previous frame
		polytope<T> P(IDSYS);
		optim_all(P,v);
		for (std::size_t X = 0; X < size_x; X++){
			traj(frame,X) = current(X);
			T dat = (1 - traj(frame,X))/TCA;           //excitation u 1
			T ddt = (0 - traj(frame,X))/TCD;           // excitation u 0
			ub(frame,X) = traj(frame,X) + DT*dat;
			lb(frame,X) = traj(frame,X) + DT*ddt;
			if (ub(frame,X) > 1) ub(frame,X) = 1; //set ubs greater than 1 to 1
			if (lb(frame,X) < 0) lb(frame,X) = 0; //set lbs less than 0 to 0
		}
	};

	void writeToCSVfile(const std::string& Aname, Eigen::MatrixXd matrix)
	{
			Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
	    ofstream file(Aname);
	    file << matrix.format(CSVFormat);
	    file.close();
	 }

	void run()
	{
		seed_trajectory();
		for (std::size_t nf = 1; nf < F; nf++){
			//cout << "frame:" << nf<< endl;
			iter(nf);
		}
	};

	void run_optim_one(std::size_t m, T v)
	{
		seed_trajectory_optim_one(m,v);
		for (std::size_t nf = 1; nf < F; nf++){
		//cout << "frame:" << nf<< endl;
			iter_optim_one(nf,m,v);
		}
	};

	void run_optim_all(T v)
	{
		seed_trajectory_optim_all(v);
		for (std::size_t nf = 1; nf < F; nf++){
			//cout << "frame:" << nf<< endl;
			iter_optim_all(nf,v);
		}
	};

	void printtraj(const std::string& fname)
	{
		cout << endl << "saving to " << fname <<endl;
		traj.print(fname);
	};
	void printlb(const std::string& fname)
	{
		cout << endl << "saving to " << fname <<endl;
		lb.print(fname);
	};
	void printub(const std::string& fname)
	{
		cout << endl << "saving to " << fname <<endl;
		ub.print(fname);
	};
	void print_traj()
	{
		cout << "traj:[" << traj.size_y << "," << traj.size_x << "]" << endl;
		traj.print();
	};
	void print_traj(std::size_t Y)
	{
		cout << "traj:frame:" << "Row(" << Y << ")@[" << traj.size_y << "," << traj.size_x << "]" << endl;
		traj.row(Y).print();
	}
	void print_lb()
	{
		cout << "lb:" << endl;
		lb.print();
	};
	void print_IDSYS(mat<T>& IDSYS)
	{
		cout << "IDSYS:[" << IDSYS.size_y << "," << IDSYS.size_x << "]" << endl;
		IDSYS.print();
	};
	void print_ub()
	{
		cout << "ub:" << endl;
		ub.print();
	};
	void print_insides()
	{
		cout << "insides:" << endl;
		insides.print();
	};
};
#endif
