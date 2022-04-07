#ifndef REDUND_H_
#define REDUND_H_
#include "base.h"
#include "util.h"
#include "kvp.h"
#include "vec.h"
#include "mat.h"
#include "polytope.h"
#include "optim.h"
//wrapper class that repeatedly calls optim
//to determine if inequality in H rep is redundant
using namespace std;
template <typename T>
class redund
{
public:
	mat<T> ActualH;
	mat<T> ActualA;
	vec<T> Actualb;
	redund(const std::string& Hin)
	{
		ActualH = mat<T>(Hin);
		ActualA = ActualH.trunc(0,1,ActualH.size_y,ActualH.size_x-1);
		Actualb = ActualH.col(0);
	};
	redund(const mat<T>& Hin)
	{
		ActualH = mat<T>(Hin);
		ActualA = ActualH.trunc(0,1,ActualH.size_y,ActualH.size_x-1);
		Actualb = ActualH.col(0);
	};
	mat<T> run()
	{
		mat<T> H= ActualH;
		mat<T> A= ActualA;
		vec<T> b= Actualb;
		vec<T> c = A.row(0);
		vec<T> tempH = H.row(0);
		T v = b(0);
		vec<bool> keepers(ActualH.size_y);
		for (std::size_t Y = 0; Y < ActualA.size_y; Y++){
			tempH=H.row(0);
			H = H.trunc(1,0,H.size_y-1,H.size_x);
			c = A.row(0);
			v = b(0);
			A = A.trunc(1,0,A.size_y-1,A.size_x);
			b = b.trunc(1,b.size-1);
			optim<T> o(H,c);
			vec<T> out = o.simplex();
			if (abs(out(out.size-1)) > v){
				H=H.concat(tempH);
				A=A.concat(c);
				b=b.append(v);
			}
		}
		return H;
	};
	polytope<T> mkp()
	{
		mat<T> Hrep = run();
		polytope<T> pol(Hrep);
		return pol;
	};
};
#endif
