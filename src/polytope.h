#ifndef POLYTOPE_H_
#define POLYTOPE_H_
#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
//polytope H-rep of
//half spaces in form A x <= b
using namespace std;
template <typename T>
class polytope
{
public:
	std::size_t size_x{},size_y{};
	mat<T> A{};
	vec<T> b{};
	mat<T> AT{};
	mat<T> matb{};
	mat<T> auxilliar_points{};
	vec<T> lambdas{};
	mat<T> HT{};
	mat<T> H{};
	vec<T> hom;
	//default constructor
	polytope()
	{
		A = mat<T>(1,1);
		b = vec<T>(1);
		AT = A.transpose();
		matb = to_mat(b);
		HT = matb.concat(-AT);
		H = HT.transpose();
		size_x = A.size_x;
		size_y = A.size_y;
		auxilliar_points = mat<T>(size_y,size_x);
		lambdas = vec<T>(size_y);
		hom = vec<T>(1); hom(0)=1;
	}
	//copy constructor
	polytope(const polytope& p)
	{
		A = p.A;
		b = p.b;
		AT = p.AT;
		matb = p.matb;
		HT = p.HT;
		H = p.H;
		size_x = A.size_x;
		size_y = A.size_y;
		auxilliar_points = p.auxilliar_points;
		lambdas = p.lambdas;
		hom = vec<T>(1); hom(0)=1;
	}
	// copy assignment constructor
	polytope& operator=(const polytope& p)
	{
		A = p.A;
		b = p.b;
		AT = p.AT;
		matb = p.matb;
		HT = p.HT;
		H = p.H;
		size_x = A.size_x;
		size_y = A.size_y;
		auxilliar_points = p.auxilliar_points;
		lambdas = p.lambdas;
		hom = vec<T>(1); hom(0)=1;
		return *this;
	}
	//A mat and b vec in constructor
	polytope(const mat<T>& Ain,const vec<T>& bin)
	{
		A = Ain;
		b = bin;
		AT = A.transpose();
		matb = to_mat(b);
		HT = matb.concat(-AT);
		H = HT.transpose();
		size_x = A.size_x;
		size_y = A.size_y;
		auxilliar_points = mat<T>(size_y,size_x);
		lambdas = vec<T>(size_y);
		hom = vec<T>(1); hom(0)=1;
	};
	//H mat constructor form [b | -A] >=0
	polytope(const mat<T>& Hin)
	{
		H = Hin;
		HT = H.transpose();
		b = H.col(0);
		A = -H.trunc(0,1,H.size_y,H.size_x-1);
		matb = to_mat(b);
		size_x = A.size_x;
		size_y = A.size_y;
		auxilliar_points = mat<T>(size_y,size_x);
		lambdas = vec<T>(size_y);
		hom = vec<T>(1); hom(0)=1;
	};
	//A mat and b vec finput constructor
	polytope(const std::string& Afile, const std::string& bfile)
	{
		A = mat<T>(Afile);
		b = vec<T>(bfile);
		AT = A.transpose();
		matb = to_mat(b);
		HT = matb.concat(-AT);
		H = HT.transpose();
		size_x = A.size_x;
		size_y = A.size_y;
		auxilliar_points = mat<T>(size_y,size_x);
		lambdas = vec<T>(size_y);
		hom = vec<T>(1); hom(0)=1;
	};
	//H finput constructor form [b | -A] >=0
	polytope(const std::string& Hfile)
	{
		H = mat<T>(Hfile);
		b = H.col(0);
		A = -H.trunc(0,1,H.size_y,H.size_x-1);
		size_x = A.size_x;
		size_y = A.size_y;
		auxilliar_points = mat<T>(size_y,size_x);
		lambdas = vec<T>(size_y);
		hom = vec<T>(1); hom(0)=1;
	};
	//destructor
	~polytope(){};
	// see if p satifies each inequality
	// fukuda says this is a bad way to do it, but is there a better way?
	bool check_inside(const vec<T>& p)
	{
		std::size_t Y=0;
		vec<T> this_row(size_y);
		for (Y = 0; Y < size_y; Y++){
			this_row = H.row(Y);
			T val((hom.concat(p))%this_row);
			if (val < -TOL){return false;}
		}
		return true;
	};
	// find auxilliars, these are points in the halfspace H_i
	void find_auxilliar_points()
	{
			std::size_t X,Y;
			auxilliar_points.set_all(0);
			for (Y = 0; Y < A.size_y; Y++){
				X = 0;
				while (abs(A(Y,X)) < TOL){
					X++;
				}
				auxilliar_points(Y,X) =b(Y)/A(Y,X);
			}
	};
	//project the auxilliars from the given point along dirvec
	//(euclidean norms to the half space along the dirvec)
	//parallel or in the halfspace returns NAN
	void compute_lambdas(const vec<T>& point,const vec<T>& direction)
	{
		std::size_t Y;
		T n,d;
		for (Y=0;Y<size_y; Y++){
			d = (direction % A.row(Y));
			if (abs(d) < TOL){
				lambdas(Y) = NAN;
			}else{
				n = (auxilliar_points.row(Y) - point) % A.row(Y);
				lambdas(Y) = n/d;
			}
		}
	};
	//+--------+
	//|PRINTERS|
	//+--------+

	void print_A()
	{
		cout << "A:[" << A.size_y << "," << A.size_x << "]" << endl;
		A.print();
	};
	void print_b()
	{
		cout << "b:[" << b.size << "]" << endl;
		b.print();
	};
	void print_aux()
	{
		cout << "Aux:[" << auxilliar_points.size_y << "," << auxilliar_points.size_x << "]" << endl;
		auxilliar_points.print();
	};
	void print_lambdas()
	{
		cout << "lam:[" << lambdas.size << "]" << endl;
		lambdas.print();
	};
	void print_H()
	{
		cout << "H:[" << H.size_y << "," << H.size_x << "]" << endl;
		H.print();
	};
	void dump()
	{
		cout << "dumping polytope:" << endl;
		print_H();
		print_A();
		print_b();
		print_aux();
		print_lambdas();
	};
#ifdef DEBINFO
	#endif
};
#endif
