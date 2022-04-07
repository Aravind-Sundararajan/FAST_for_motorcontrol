#ifndef OPTIM_H_
#define OPTIM_H_
#include "base.h"
#include "util.h"
#include "kvp.h"
#include "vec.h"
#include "mat.h"
#include "polytope.h"
using namespace std;
//optimizer class that performs simplex algorithm in input linear program.
//can solve mixed constraint problems Ax<=b or Ax >= b
//form [b|-A]>=0
//decides to add artificial and slack variables based on inequality direction
//phase 1 problem: minimize the artificial variables
//phase 2 problem: optimize the phase 1 solution using original objective
template <typename T>
class optim
{
public:
	std::size_t size_x, size_y, artc, state, iter;
	mat<T> H,A,constr,tableau;
	vec<T> b,c,c_old,c_og,Xsol;
	vec<std::size_t> ai;
	bool twophase = false;
	T v;
	//constructor using A mat and b vec, form [b|-A] >= 0
	optim(const mat<T>& Ain,const vec<T>& bin,const vec<T>& cin)
	{
		A = Ain;
		b = bin;
		c = cin;
		T cc = c(0);
		c = c.trunc(1,c.size-1);
		v=0;
		ai =vec<std::size_t>(1);
		ai(0)=0;
		c_old = c;
		c_og = c;
		size_x = Ain.size_x;
		size_y = Ain.size_y;
		state = 0;
		iter = 0;
		A = A.transpose();
		for (std::size_t Y = 0; Y < size_y; Y++){
			if ((b(Y) < 0) || (b(Y) == -0)){
				vec<T> temp(size_y);
				temp(Y) = -1;
				A=A.concat(temp);
				c=c.append(0);
			}
			else{
				if (b(Y) >= 0){
					vec<T> temp(size_y);
					temp(Y) = 1;
					A=A.concat(temp);
					c=c.append(0);
				}
			}
		}
		for (std::size_t Y = 0; Y < size_y; Y++){
			if ((b(Y) < 0) || (b(Y) == -0)){
				twophase=true;
				artc++;
				b(Y) = -b(Y);
				vec<T> temp(size_y);
				temp(Y) = 1;
				ai=ai.append(Y);
				A=A.concat(temp);
				c=c.append(-1);
			}
		}
		vec<T> temp(size_y);
		constr = (A.concat(b)).transpose();
		tableau = to_mat(-c.append(cc)).concat(constr);
		vec<T> zj(tableau.size_x);
		tableau = tableau.concat(zj); //zj row
		tableau = tableau.concat(tableau.row(tableau.size_y-2)); // Cj- zj row
	};
	// constructor using H mat
	optim(const mat<T>& Hin,const vec<T>& cin)
	{
		H = Hin;
		A = -H.trunc(0,1,H.size_y,H.size_x-1);
		b = H.col(0);
		c = cin;
		T cc = c(0);
		c = c.trunc(1,c.size-1);
		ai =vec<std::size_t>(1);
		ai(0)=0;
		v=0;
		c_old = c;
		c_og = c;
		size_x = A.size_x;
		size_y = A.size_y;
		state = 0;
		iter = 0;
		A = A.transpose();
		for (std::size_t Y = 0; Y < size_y; Y++){
			if ((b(Y) < 0) || (b(Y) == -0)){
				vec<T> temp(size_y);
				temp(Y) = -1;
				A=A.concat(temp);
				c=c.append(0);
			}
			else{
				if (b(Y) >= 0){
					vec<T> temp(size_y);
					temp(Y) = 1;
					A=A.concat(temp);
					c=c.append(0);
				}
			}
		}
		for (std::size_t Y = 0; Y < size_y; Y++){
			if ((b(Y) < 0) || (b(Y) == -0)){
				twophase=true;
				artc++;
				b(Y) = -b(Y);
				vec<T> temp(size_y);
				temp(Y) = 1;
				ai=ai.append(Y);
				A=A.concat(temp);
				c=c.append(-1);
			}
		}
		vec<T> temp(size_y);
		constr = (A.concat(b)).transpose();
		tableau = to_mat(-c.append(cc)).concat(constr);
	};
	//copy constructor
	optim(const optim<T>& o)
	{
		size_x = o.size_x;
		size_y = o.size_y;
		artc = o.artc;
		state = o.state;
		iter = o.iter;
		H = o.H;
		A = o.A;
		constr = o.constr;
		tableau = o.tableau;
		b = o.b;
		c = o.c;
		c_old = o.c_old;
		c_og = o.c_og;
		Xsol = o.Xsol;
		ai = o.ai;
		twophase = o.twophase;
		v = o.v;
	};
	//copy assignment operator
	optim<T>& operator=(const optim<T>& o)
	{
		size_x = o.size_x;
		size_y = o.size_y;
		artc = o.artc;
		state = o.state;
		iter = o.iter;
		H = o.H;
		A = o.A;
		constr = o.constr;
		tableau = o.tableau;
		b = o.b;
		c = o.c;
		c_old = o.c_old;
		c_og = o.c_og;
		Xsol = o.Xsol;
		ai = o.ai;
		twophase = o.twophase;
		v = o.v;
		return *this;
	};
	//destructor
	~optim(){};
	//debinfo method for printing the current simplex tableau
	#ifdef DEBINFO
	void print_tableau()
	{
		cout << "ITER:" << iter << endl;
		cout << "size_y:" << tableau.size_y << ",size_x:"<< tableau.size_x << endl<< setfill(' ') << setw(4) << " ";
		for (std::size_t X =0; X < tableau.size_x-1; X++) {
			cout << setprecision(3) << setfill(' ') << setw(7) << X << "   ";
		}
		cout << endl;
		cout << setfill(' ') << setw(4) <<"-+";
		for (std::size_t X =0; X < tableau.size_x-1; X++) {
			cout << setprecision(3) << setfill('-') << setw(7) << '-' << "---";
		}
		cout << '+'<< endl << setfill(' ') << setw(4) << "C|";
		for (std::size_t X =0; X < tableau.size_x-1; X++) {
			cout << setprecision(3) << setfill(' ') << setw(7) << tableau(0,X) << "   ";
		}
		cout <<'|' << setprecision(3) << setfill(' ') << setw(7) << tableau(0,tableau.size_x-1) << endl<< setfill(' ') << setw(4) <<"-+";
		for (std::size_t X =0; X < tableau.size_x-1; X++) {
			cout << setprecision(3) << setfill('-') << setw(7) << '-' << "---";
		}
		cout <<'+' << setprecision(3) << setfill('-') << setw(7) << '-' << endl;
		for (std::size_t Y =1; Y < tableau.size_y-2; Y++) {
			cout<< setfill(' ') << setw(3) << Y << "|";
			for (std::size_t X =0; X < tableau.size_x-1; X++) {
				cout << setprecision(3) << setfill(' ') << setw(7) << tableau(Y,X) << "   ";
			}
			cout <<'|' << setprecision(3) << setfill(' ') << setw(7) << tableau(Y,tableau.size_x-1) << endl;
		}
		cout << setfill(' ') << setw(4) <<"-+";
		for (std::size_t X =0; X < tableau.size_x-1; X++) {
			cout << setprecision(3) << setfill('-') << setw(7) << '-' << "---";
		}
		cout << endl;
		for (std::size_t Y =tableau.size_y-2; Y < tableau.size_y; Y++) {
			cout<< setfill(' ') << setw(3) << Y << "|";
			for (std::size_t X =0; X < tableau.size_x-1; X++) {
				cout << setprecision(3) << setfill(' ') << setw(7) << tableau(Y,X) << "   ";
			}
			cout <<'|' << setprecision(3) << setfill(' ') << setw(7) << tableau(Y,tableau.size_x-1) << endl;
		}
		cout << endl;
	};
	#endif
	//identify pivot column
	//this is the column with the objective row element that is the most negative
	kvp<std::size_t,T> find_pivot_col()
	{
		std::size_t cp = OOB;
		T cv = BIGM;
		for (std::size_t X =0; X < tableau.size_x-1; X++){
			if (tableau(0,X) < cv) {
				cv = tableau(0,X);
				cp = X;
			}
		}
		return kvp<std::size_t,T>(cp,cv);
	};
	//identify pivot row
	//take the elementwise quotient of the pivot column and the b col
	//the least positive value is the pivot row
	kvp<std::size_t,T> find_pivot_row(std::size_t X)
	{
		std::size_t rp = OOB;
		T rv = BIGM;
		vec<T> pivot_column = tableau.col(tableau.size_x-1) / tableau.col(X);
		for (std::size_t Y =1 ; Y < tableau.size_y; Y++){
			if (((pivot_column(Y) < rv) && (pivot_column(Y) > 0))){
				rv = pivot_column(Y);
				rp = Y;
			}
		}
		return kvp<std::size_t,T>(rp,rv);
	};
	void pivot(std::size_t Y, std::size_t X)// pivot Yth variable around xth constraint
	{
		tableau = tableau.divrow(Y,tableau(Y,X));
		#ifdef DEBINFO
		cout << "Pivoting x_"<< X << " around constraint " << Y << "." << endl;
		print_tableau();
		cout << "Row ops, turning basic variable into nonbasic." << endl;
		#endif
		for (std::size_t YY =0; YY < tableau.size_y; YY++){
			if (YY != Y){
				#ifdef DEBINFO
				cout << "R("<<YY<< ")->" << "R("<<YY<<") +"<<-tableau(YY,X)<< "" <<"*R("<<Y<<")" << endl;
				#endif
				tableau = tableau.rowop(YY,Y,-tableau(YY,X));
			}
		}
		#ifdef DEBINFO
		print_tableau();
		#endif
	};
	std::size_t iter_simplex()// Returns 0 if OK, -1 if INFEASIBLE
	{
		bool oc = true;
		for (std::size_t X =0; X < tableau.size_x; X++){
			if (tableau(0,X) < 0) oc = false;
		}
		if (oc)	return 0;
		iter++;
		kvp<std::size_t,T> colp = find_pivot_col();
		kvp<std::size_t,T> rowp = find_pivot_row(colp.first);
		if ((rowp.second != std::numeric_limits<T>::max())  && (colp.second != std::numeric_limits<T>::max())){
			if ((rowp.first != OOB) && (colp.first != OOB)){
				pivot(rowp.first,colp.first);
				for (std::size_t X =0; X < tableau.size_x-1; X++){
					if (tableau(0,X) < 0){return 3;}//not done.
				}
				return 0; // done.
			}
			else{return 1;} // infeasible problem!
		}else{return 2;} //the problem is unbounded
		return OOB;
	};
	// Runs the simplex algorithm to optimise the LP.
	vec<T> simplex()//2  if unbounded, 1 if infeasible,3 if not done,OOB if done.
	{
		#ifdef DEBINFO
		cout<< "PROBLEM:" << endl;
		print_tableau();
		#endif
		if (!twophase){
			//rewrite as nonbasics
			// for (std::size_t X = 0; X < A.size_x; X++){
			// 	if (tableau(0,X) !=0){
			// 		std::size_t Y =0;
			// 		while(true){
			// 			Y++;
			// 			if (tableau(Y,X) !=0){
			// 				tableau = tableau.rowop(0,Y,-tableau(0,X)/tableau(Y,X));
			// 				break;
			// 			}
			// 		}
			// 	}
			// }
			#ifdef DEBINFO
			cout<< "rewriting as nonbasics:" << endl;
			print_tableau();
			#endif
			while ((state = iter_simplex()) == 3);
			if (state == 1){
				vec<T> ret = tableau.row(0);
				ret.set_all(-1); // infeasible
				return ret;
			}
			else if (state == 2){
				vec<T> ret = tableau.row(0);
				ret.set_all(-2); //unbounded
				return ret;
			}
			vec<T> ret = (tableau.row(0));
			#ifdef DEBINFO
			cout << "optimal solution:" << endl;
			ret.print();
			cout << "optimal value:" << tableau(0,tableau.size_x-1)<<  endl;
			#endif
			return ret;
		}else{
			#ifdef DEBINFO
			cout << "the problem is twophase with "<< artc << " artificial vars." << endl;
			#endif
			c_old = tableau.row(0);
			for (std::size_t X = 0; X < tableau.size_x-artc-1; X++){
				tableau(0,X) = 0;
			}
			tableau(0,size_x-1) = 0;
			#ifdef DEBINFO
			cout<< "MODIFIED PHASE 1 PROBLEM:" << endl;
			print_tableau();
			#endif
			//rewrite as nonbasics
			for (std::size_t Y = 1; Y < ai.size; Y++){
				tableau = tableau.rowop(0,ai(Y)+1,-1);
			}
			#ifdef DEBINFO
			ai.print();
			cout<< "REWRITE AS NONBASIC:" << endl;
			print_tableau();
			#endif
			while (iter_simplex() == 3); //PHASE 1
			if (state == 1){
				vec<T> ret = tableau.row(0);
				ret.set_all(-1); // infeasible
				return ret;
			}
			else if (state == 2){
				vec<T> ret = tableau.row(0);
				ret.set_all(-2); //unbounded
				return ret;
			}
			vec<T> ret = tableau.row(0);
			#ifdef DEBINFO
			cout << "phase 1 solution:" << endl;
			ret.print();
			cout << "phase 1 optimal:" << tableau(0,tableau.size_x-1)<<  endl;
			#endif
			iter = 0;
			for (std::size_t X = 0; X < tableau.size_x-artc; X++){
				tableau(0,X) = c_old(X);
			}
			mat<T> b_old = tableau.trunc(0,tableau.size_x-1,tableau.size_y,1);
			mat<T> A_old = tableau.trunc(0,0,tableau.size_y,tableau.size_x-artc-1);
			tableau = (A_old.transpose()).concat(b_old.transpose()).transpose();
			#ifdef DEBINFO
			cout<< "PHASE 2 PROBLEM:" << endl;
			print_tableau();
			#endif
			//rewrite as nonbasics
			for (std::size_t X = 1; X < ai.size; X++){
				if (tableau(0,ai(X))!=0){
					T mod = -dbl::max();
					std::size_t modrow;
					for (std::size_t Y = 1; Y < tableau.size_y; Y++){
						if ((tableau(Y,ai(X)) > mod) && (tableau(Y,ai(X)) != 0)){
							mod = tableau(Y,ai(X));
							modrow = Y;
						}
					}
					tableau = tableau.rowop(0,modrow,-tableau(0,ai(X))/mod);
				}
			}
			#ifdef DEBINFO
			cout << "REWRITE AS NONBASIC:" << endl;
			print_tableau();
			#endif
			while (iter_simplex() == 3); //PHASE 2
			if (state == 1){
				vec<T> ret = tableau.row(0);
				ret.set_all(-1); // infeasible
				return ret;
			}
			else if (state == 2){
				vec<T> ret = tableau.row(0);
				ret.set_all(-2); //unbounded
				return ret;
			}
			ret= tableau.row(0);
			#ifdef DEBINFO
			cout << "phase 2 solution::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
			ret.print();
			cout << "phase 2 optimal:" << tableau(0,tableau.size_x-1)<<  endl;
			#endif
			return ret;
		}
	};
};
#endif
