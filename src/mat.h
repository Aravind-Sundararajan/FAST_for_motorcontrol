#ifndef MAT_H_
#define MAT_H_
#include "base.h"
#include "util.h"
#include "kvp.h"
#include "vec.h"
//container template class that holds a matrix as a 1D array
//convenient operators, arbitrary function application on elem
//rowops, printers, etc...
using namespace std;

template <typename T>
class mat
{
private:

	//returns the determinant of a matrix, wrapper function det() calls this
	T determinant(mat<T>& temp, std::size_t n)
	{
		T det=0;
		std::size_t p, h, k, i, j;
		if(n==1) {
			return get(0,0);
		} else if(n==2) {
			det=(get(0,0) * get(1,1) - get(0,1) * get(1,0));
			return det;
		} else {
			for(p=0;p<n;p++) {
				h = 0;
				k = 0;
				for(i=1;i<n;i++) {
					for( j=0;j<n;j++) {
						if(j==p) {
							continue;
						}
						temp(h,k) = get(i,j);
						k++;
						if(k==n-1) {
							h++;
							k = 0;
						}
					}
				}
				det = det + get(0,p) * pow(-1,p)*determinant(temp,n-1);
			}
			return det;
		}
	};
public:
	T *A{};
	std::size_t size_x{};
	std::size_t size_y{};
	//default constructor
	mat()
	{
		size_x = (std::size_t)1;
		size_y = (std::size_t)1;
		A = new T [size_x*size_y]{};
		*(A) = (T)0;
	};
	//size constructor, makes a zero matrix of size sy * sx
	mat(std::size_t sy,std::size_t sx)
	{
		size_x = sx;
		size_y = sy;
		A = new T [size_x*size_y]{};
	};
	//finput constructor
	mat(const std::string& Aname)
	{
		size_x = nCol(Aname);
		size_y = nRow(Aname);
		A = new T [size_x*size_y]{};
		set_all(Aname);
	};
	//copy constructor
	mat(const mat<T>& v)
	{
		size_x = v.size_x;size_y = v.size_y;
		A = new T [size_x*size_y]{};
		for (std::size_t Y =0; Y< size_y; Y++){
			for (std::size_t X =0; X< size_x; X++){
				A[Y+size_y*X] = v(Y,X);
			}
		}
	};
	//copy assignment
	mat<T>& operator=(const mat<T>& v){
		size_x = v.size_x;size_y = v.size_y;
		delete[] A;
		A = new T [size_x*size_y]{};
		for (std::size_t Y =0; Y< size_y; Y++){
			for (std::size_t X =0; X< size_x; X++){
				A[Y+size_y*X] = v(Y,X);
			}
		}
		return *this;
	}
	//destructor
	~mat()
	{
		delete[] A;
	};
	//dynamic casting
	template <typename U>
	mat<T>& operator=(const mat<U>& v)
	{
		delete[] A;
		size_x = v.size_x;size_y = v.size_y;
		A = new T [size_x*size_y]{};
		for (std::size_t Y =0; Y< size_y; Y++){
			for (std::size_t X =0; X< size_x; X++){
				A[Y+size_y*X] = (T) v(Y,X);
			}
		}
		return *this;
	}
	//template cast, we're in C++ but we get MATLAB style casting
	//(this works if the cast is defined)
	// mat<int> i(size_y,size_x);
	// mat<char> c(i);
	template <typename U>
	mat(const mat<U>& v)
	{
		size_x = v.size_x;
		size_y = v.size_y;
		A = new T [size_x*size_y]{};
		for (std::size_t Y =0; Y< size_y; Y++){
			for (std::size_t X =0; X< size_x; X++){
				A[Y+size_y*X] = (T) v(Y,X);
			}
		}
	}
	//returns value at y,x
	T& operator()(std::size_t Y,std::size_t X)
	{
		return A[Y+size_y*X];
	}
	//returns const value at y,x
	T& operator()(std::size_t Y,std::size_t X) const
	{
		return A[Y+size_y*X];
	}
	//funny looking set operator m(y,x,v);
	//kinda looks like a tuple
	void operator()(std::size_t Y,std::size_t X, T v)
	{
		A[Y+size_y*X] = v;
	}
	//returns value at y,x, same as the get operator
	T& get(std::size_t Y,std::size_t X)
	{
		return A[Y+size_y*X];
	}
	//returns const value at y,x, same as the get operator
	T& get(std::size_t Y,std::size_t X) const
	{
		return A[Y+size_y*X];
	}
	//set value at y,x, same as the set operator
	void set(std::size_t Y,std::size_t X,T v)
	{
		A[Y+size_y*X] = v;
	}
	//set every index in matrix to v
	void set_all(T v)
	{
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				set(Y,X, v);
			}
		}
	}
	//copy values by index of mat
	void set_all(const mat<T>& v)
	{
		std::size_t szx, szy;
		if (v.size_x < size_x){
			szx = v.size_x;
		}else {
			szx = size_x;
		}
		if (v.size_y < size_y){
			szy = v.size_y;
		}else {
			szy = size_y;
		}
		for (std::size_t Y = 0; Y < szy; Y++){
			for (std::size_t X = 0; X < szx; X++){
				set(Y,X, v(Y,X));
			}
		}
	}

	mat<T> diag(const vec<T>& v)
	{
		mat<T> o(v.size,v.size);
		for ( std::size_t Y = 0; Y < v.size; Y++ ){
			o(Y,Y) = v(Y);
		}
		return o;
	};

	//set all from input file
	void set_all(const std::string& Aname)
	{
		string line,str;
		T f;
		ifstream Afile(Aname);
		std::size_t Y = 0;
		while(getline(Afile,line) && !line.empty()) {
			istringstream issline(line);
			std::size_t X = 0;
			while(getline(issline,str,',')){
				istringstream ss(str);
				ss >> f;
				A[Y+size_y*X] = f;
				X++;
				ss.clear();
			}
			issline.clear();
			Y++;
		}
		Afile.close();
	}
	//this does nothing
	mat<T> operator+ ()
	{
		mat<T> u(*this);
		return u;
	};
	//matrix negation
	mat<T> operator- ()
	{
		mat<T> u(*this);
		u*=(T)-1;
		return u;
	};
	//matrix addition with scalar  (element-wise)
	mat<T> operator+ (T v)
	{
		mat<T> u(*this);
		u+=v;
		return u;
	};
	//matrix subtraction with scalar  (element-wise)
	mat<T> operator- (T v)
	{
		mat<T> u(*this);
		u-=v;
		return u;
	};
	//matrix-matrix additon
	mat<T> operator+ (const mat<T>& v)
	{
		mat<T> u(*this);
		u+=v;
		return u;
	};
	//matrix-matrix subtraction
	mat<T> operator- (const mat<T>& v)
	{
		mat<T> u(*this);
		u-=v;
		return u;
	};
	//matrix scalar multiplication (element-wise)
	mat<T> operator* (T v)
	{
		mat<T> u(*this);
		u*=v;
		return u;
	};
	//matrix scalar division (element-wise)
	mat<T> operator/ (T v)
	{
		mat<T> u(*this);
		u/=v;
		return u;
	};
	//matrix-matrix multiplication  (element-wise)
	mat<T> operator* (const mat<T>& v)
	{
		mat<T> u(*this);
		u*=v;
		return u;
	};
	//matrix-matrix division  (element-wise)
	mat<T> operator/ (const mat<T>& v)
	{
		mat<T> u(*this);
		u/=v;
		return u;
	};
	//matmul A_ay,ax % B_by,bx = C_ay,bx
	mat<T> operator% (const mat<T>& v)
	{
		std::size_t i,j,k,r1 = size_y,	c1 = size_x,	c2 = v.size_x;
		mat<T> out(r1,c2);//r2 = v.size_y
		for ( i = 0; i < r1; i++){
			for (j = 0; j < c2; j++){
				for (k = 0; k < c1; k++){
					out(i,j)+= (*this)(i,k) * v(k,j);
				}
			}
		}
		return out;
	};
	//multiply row y with scalar
	mat<T> multrow (std::size_t Y,T v)
	{
		mat<T> out(*this);
		for (std::size_t X = 0; X < size_x; X++){
			out(Y,X)=out(Y,X) * v;
		}
		return out;
	}
	//divide row Y with scalar
	mat<T> divrow (std::size_t Y,T v)
	{
		mat<T> out(*this);
		for (std::size_t X = 0; X < size_x; X++){
			out(Y,X)=out(Y,X) / v;
		}
		return out;
	}
	// add row Y2 to row Y
	// R(Y) -> R(Y) + R(Y2)
	mat<T> rowop (std::size_t Y,std::size_t Y2)
	{
		mat<T> out(*this);
		for (std::size_t X = 0; X < size_x; X++){
			out(Y,X)=out(Y,X) + out(Y2,X);
		}
		return out;
	};
	// add row Y2 multiplied by scalar to row Y
	// R(Y) -> R(Y) + v * R(Y2)
	mat<T> rowop (std::size_t Y,std::size_t Y2, T v)
	{
		mat<T> out(*this);
		for (std::size_t X = 0; X < size_x; X++){
			out(Y,X)=out(Y,X) + v*out(Y2,X);
		}
		return out;
	};
	// add row Y2 to row Y
	// R(Y) -> R(Y) + v
	mat<T> rowop (std::size_t Y,const vec<T>& v)
	{ // add v to row Y: R(Y) = R(Y) + v(Y2)
		mat<T> out(*this);
		for (std::size_t X = 0; X < size_x; X++){
			out(Y,X)=out(Y,X) + v(X);
		}
		return out;
	};
	//swap row Y with row Y2
	//temp  -> R(Y)
	//R(Y)  -> R(Y2)
	//R(Y2) -> temp
	mat<T> rowswap (std::size_t Y,std::size_t Y2)
	{
		mat<T> out(*this);
		vec<T> temp(out.row(Y));
		for (std::size_t X = 0; X < size_x; X++){
			out(Y,X)= out(Y2,X);
			out(Y2,X) = temp(X);
		}
		return out;
	};
	//assignment scalar addition operator (elementwise)
	mat<T>& operator+=(T v)
	{
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				A[Y+size_y*X] = A[Y+size_y*X] + v;
			}
		}
		return *this;
	};
	//assignment scalar subtraction (elementwise)
	mat<T>& operator-=(T v)
	{
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				A[Y+size_y*X] = A[Y+size_y*X] - v;
			}
		}
		return *this;
	};
	//assignment matrix addition (elementwise)
	mat<T>& operator+=(const mat<T>& v)
	{
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				A[Y+size_y*X] = A[Y+size_y*X] + v(Y,X);
			}
		}
		return *this;
	};
	//assignment matrix subtraction (elementwise)
	mat<T>& operator-=(const mat<T>& v)
	{
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				A[Y+size_y*X] = A[Y+size_y*X] - v(Y,X);
			}
		}
		return *this;
	};
	//assignment scalar multiplication (elementwise)
	mat<T>& operator*=(T v)
	{
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				A[Y+size_y*X] = A[Y+size_y*X] * v;
			}
		}
		return *this;
	};
	//assignment scalar division (elementwise)
	mat<T>& operator/=(T v)
	{
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				A[Y+size_y*X] = A[Y+size_y*X] / v;
			}
		}
		return *this;
	};
	//assignment matrix multiplication (elementwise)
	mat<T>& operator*=(const mat<T>& v)
	{
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				A[Y+size_y*X] = A[Y+size_y*X] * v(Y,X);
			}
		}
		return *this;
	};
	//assignment matrix division (elementwise)
	mat<T>& operator/=(const mat<T>& v)
	{
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				A[Y+size_y*X] = A[Y+size_y*X] / v(Y,X);
			}
		}
		return *this;
	};
	//mat fun, do f() to every element
	mat<T> fun(T f(T)) // pass function pointer
	{
		mat<T> out(size_y,size_x);
		std::size_t X,Y;
		for (Y =0; Y< size_y; Y++){
			for (X =0; X< size_x; X++){
				out(Y,X)=f((*this)(Y,X));
			}
		}
		return out;
	};
	//sum rows or columns and return vec
	vec<T> sum(bool dim)//0,1
	{
		vec<T> out;
		std::size_t X,Y;
		if (dim){
			out = vec<T>(size_x);
			for (Y =0; Y< size_y; Y++){
				for (X =0; X< size_x; X++){
					out(X)+=(*this)(Y,X);
				}
			}
		}
		else{
			out = vec<T>(size_y);
			for (Y =0; Y< size_y; Y++){
				for (X =0; X< size_x; X++){
					out(Y)+=(*this)(Y,X);
				}
			}
		}
		return out;
	};
	//returns the matrix adjoint or adjunct
	mat<T> adjoint()
	{
		std::size_t xi=0,yi=0,X=0,Y=0,xid=0,yid=0, r1 = size_y,c1 = size_x,r2=r1-1,c2=c1-1;
		mat<T> adj(r1,c1);
		adj.set_all(0);
		mat<T> cof(r2,c2);
		for ( yi = 0; yi < r1; yi++){
			for ( xi = 0; xi < c1; xi++){
				cof.set_all(0);
				yid=0;
				for ( Y = 0; Y < r1; Y++){
					xid=0;
					if (Y!=yi){
						for ( X = 0; X < r1; X++){
							if (X!=xi){
								cof(yid,xid) = A[Y+size_y*X];
								xid++;
							}
						}
						yid++;
					}
				}
				adj(yi,xi) = pow(-1,xi)*pow(-1,yi)*cof.det();
			}
		}
		return adj.transpose();
	};
	//returns matrix inverse  inverse = adjoint/ determinant
	mat<T> inverse()
	{
		T d =det();
		mat<T> adj(adjoint());
		mat<T> inv(adj/d);
		return inv;
	};

	T frobenius()
	{//computes the frobenius norm of the matrix (sum squares all elem)
		double value;
		value = 0.0;
		for ( std::size_t Y = 0; Y < size_y; Y++ ){
			for (std::size_t X = 0; X < size_x; X++ ){
				value = value + pow ( get(Y,X), 2 );
			}
		}
		value = sqrt (value);
		return value;
	}
	T is_eigen_right(const mat<T>& x, const vec<T>& lambda)
	{
		//(*this) should be [N x N]
		// x should be [N x K]
		//lambda should be [1 x K]
		std::size_t i,j,k,l;
		std::size_t n = x.size_x;
		k = lambda.size;
		mat<T> c(n,k);
		for ( j = 0; j < k; j++ ){
			for ( i = 0; i < n; i++ ){
				c(i,j) = 0.0;
				for ( l = 0; l < n; l++ ){
					c(i,j) = c(i,j) + get(i,l) * x(i,j);
				}
			}
		}
		for ( j = 0; j < k; j++ ){
			for ( i = 0; i < n; i++ ){
				c(i,j) = c(i,j) - lambda(j) * x(i,j);
			}
		}
		T ef = frobenius ( n, k, c );
		return ef;
	};
	void identity()
	{
		(*this).set_all(0);
		for ( std::size_t Y = 0; Y < size_y; Y++ ){
			(*this)(Y,Y) = 1;
		}
	};
	vec<T> get_diag()
	{
		vec<T> v(size_y);
		for ( std::size_t Y = 0; Y < size_y; Y++ ){
			v(Y) = (*this)(Y,Y);
		}
		return v;
	}

	mat<T> jacobi_eigen()
	{
		std::size_t n =size_y;
		vec<T> bw(size_y);
    T c;
    T g;
    T gapq;
    T h;
    std::size_t i;
    std::size_t j;
    std::size_t k;
    std::size_t l;
    std::size_t m;
    std::size_t p;
    std::size_t q;
    T s;
    T t;
    T tau;
    T term;
    T termp;
    T termq;
    T theta;
    T thresh;
    T w;
    vec<T> zw(n);
		mat<T> v(n,n);
		vec<T> d(n);
		mat<T> a(*this);
		v.identity();
		d = (*this).get_diag();
		for (i = 0; i < n; i++ ){
			bw(i) = d(i);
			zw(i) = 0.0;
		}
		std::size_t it_num = 0;
		std::size_t rot_num = 0;
		while ( it_num < MAXITER ){
      it_num = it_num + 1;
      //
      //  The convergence threshold is based on the size of the elements in
      //  the strict upper triangle of the matrix.
      //
      thresh = 0.0;
      for ( j = 0; j < n; j++ ){
        for ( i = 0; i < j; i++ ){
          thresh = thresh + a(i,j) * a(i,j);
        }
      }
      thresh = sqrt ( thresh ) / ( double ) ( 4 * n );
      if ( thresh == 0.0 ){
        break;
      }
      for ( p = 0; p < n; p++ ){
        for ( q = p + 1; q < n; q++ ){
          gapq = 10.0 * fabs ( a(p,q) );
          termp = gapq + fabs ( d(p) );
          termq = gapq + fabs ( d(q) );
          //
          //  Annihilate tiny offdiagonal elements.
          //
          if ( 4 < it_num && termp == fabs ( d(p) ) && termq == fabs ( d(q) ) ){
              a(p,q) = 0.0;
            }
            //
            //  Otherwise, apply a rotation.
            //
            else if ( thresh <= fabs ( a(p,q) ) ){
              h = d(q) - d(p);
              term = fabs ( h ) + gapq;
              if ( term == fabs ( h ) ){
                t = a(p,q) / h;
              }
              else
              {
                theta = 0.5 * h / a(p,q);
                t = 1.0 / ( fabs ( theta ) + sqrt ( 1.0 + theta * theta ) );
                if ( theta < 0.0 ){
                  t = - t;
                }
              }
              c = 1.0 / sqrt ( 1.0 + t * t );
              s = t * c;
              tau = s / ( 1.0 + c );
              h = t * a(p,q);
              //
              //  Accumulate corrections to diagonal elements.
              //
              zw(p) = zw(p) - h;
              zw(q) = zw(q) + h;
              d(p) = d(p) - h;
              d(q) = d(q) + h;
              a(p,q) = 0.0;
              //
              //  Rotate, using information from the upper triangle of A only.
              //
              for ( j = 0; j < p; j++ ){
                g = a(j,p);
                h = a(j,q);
                a(j,p) = g - s * ( h + g * tau );
                a(j,q) = h + s * ( g - h * tau );
              }
              for ( j = p + 1; j < q; j++ ){
                g = a(p,j);
                h = a(j,q);
                a(p,j) = g - s * ( h + g * tau );
                a(j,q) = h + s * ( g - h * tau );
              }
              for ( j = q + 1; j < n; j++ ){
                g = a(p,j);
                h = a(q,j);
                a(p,j) = g - s * ( h + g * tau );
                a(q,j) = h + s * ( g - h * tau );
              }
              //
              //  Accumulate information in the eigenvector matrix.
              //
              for ( j = 0; j < n; j++ ){
                g = v(j,p);
                h = v(j,q);
                v(j,p) = g - s * ( h + g * tau );
                v(j,q) = h + s * ( g - h * tau );
              }
              rot_num = rot_num + 1;
            }
          }
        }
        for ( i = 0; i < n; i++ ){
          bw(i) = bw(i) + zw(i);
          d(i) = bw(i);
          zw(i) = 0.0;
        }
      }
      //
      //  Restore upper triangle of input matrix.
      //
      for ( j = 0; j < n; j++ ){
        for ( i = 0; i < j; i++ ){
          a(i,j) = a(j,i);
        }
      }
      //
      //  Ascending sort the eigenvalues and eigenvectors.
      //
      for ( k = 0; k < n - 1; k++ ){
        m = k;
        for ( l = k + 1; l < n; l++ ){
          if ( d(l) < d(m) ){
            m = l;
          }
        }
        if ( m != k ){
          t    = d(m);
          d(m) = d(k);
          d(k) = t;
          for ( i = 0; i < n; i++ ){
            w        = v(i,m);
            v(i,m) = v(i,k);
            v(i,k) = w;
          }
        }
      }
		v = v.concat(diag(d));
		return v;
	}

	//wrapper function that reurns determinant
	T det() //wrapper
	{
		mat<T> a(size_y,size_x);
		a.set_all(*this);
		T d =determinant(a,size_x);
		return d;
	}
	//returns matrix transpose
	mat<T> transpose()
	{
		mat<T> o(size_x,size_y);
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				o(X,Y) = get(Y,X);
			}
		}
		return o;
	};
	//returns row of matrix as vector
	vec<T> row(std::size_t Y) const
	{
		vec<T> out(size_x);
		for (std::size_t X =0; X< size_x; X++){
			out(X) = get(Y,X);
		}
		return out;
	}
	//returns col of matrix as vector
	vec<T> col(std::size_t X) const
	{
		vec<T> out(size_y);
		for (std::size_t Y =0; Y< size_y; Y++){
			out(Y) = get(Y,X);
		}
		return out;
	};
	//concatenates vector as row
	mat<T> concat(const vec<T>& b)  const
	{
		std::size_t X;
		std::size_t new_y = size_y+1;
		std::size_t new_x = size_x;
		mat<T> out(new_y,new_x);
		out.set_all(*this);
		for (X =0; X < size_x; X++){
			out(size_y,X) = b(X);
		}
		return out;
	};
	//concatenates rows of input matrix as rows
	mat<T> concat(const mat<T>& b) const
	{
		std::size_t Y;
		mat<T> out(size_y,size_x);
		out.set_all(*this);
		for (Y =0; Y < b.size_y; Y++){
			out = out.concat(b.row(Y));
		}
		return out;
	};
	//returns matrix without the specified row
	mat<T> cutrow(std::size_t Y)
	{
		mat<T> out(1,size_x);
		for (std::size_t YY =0; YY< size_y; YY++){
			if (YY!=Y){
				out=out.concat((*this).row(YY));
			}
		}
		out=out.trunc(1,0,size_y-1,size_x);
		return out;
	};
	//truncates matrix from starting location Y0,X0 to size szy,szx
	mat<T> trunc(std::size_t Y0,std::size_t X0,std::size_t szy,std::size_t szx)
	{
		mat<T> out(szy,szx);
		for (std::size_t Y =Y0; Y< Y0+szy; Y++){
			for (std::size_t X =X0; X< X0+szx; X++){
				out(Y-Y0,X-X0,get(Y,X));
			}
		}
		return out;
	};
	////////////////////////////////////////////
	//this following stuff is super hacky
	//and im duct-taping lrs into this stuff via extern "C"
	//so we aren't lame and shellexec()/system() calls
	//no memory leaks but,I'm not a C/C++ programmer.
	//Someone find a CS undergrad....
	///////////////////////////////////////////
	//tweaking some code from lrslib.h
	/* generate H-representation lrs from input mat<T> */
	/* construct the lrs linear problem from the input mat<T> H and objective vec<T> c*/
	void makeLP (lrs_dic *P, lrs_dat *Q,const mat<T>& M,const vec<T>& c)
	{
		long row, j;
		long m=Q->m;       /* number of inequalities      */
		//long n=Q->n;       /* hypercube has dimension n-1 */
		long d=P->d;
		lrs_mp_vector Num, Den;
		for (row=1;row<=m;row++){
			Num=lrs_alloc_mp_vector(d+1);
			Den=lrs_alloc_mp_vector(d+1);
			for(j=0;j<=d;j++){
				char num[MAXCOL];
				char den[MAXCOL];
				double t = M(row-1,j);
				rat_approx(t,num,den);
				atomp(num,Num[j]);
				atomp(den,Den[j]);
			}
			lrs_set_row_mp(P,Q,row,Num,Den,GE);
			lrs_clear_mp_vector(Num,d+1);
			lrs_clear_mp_vector(Den,d+1);
		}

		Num=lrs_alloc_mp_vector(d+1);
		Den=lrs_alloc_mp_vector(d+1);
		for(j=0;j<=d;j++){
			char numc[MAXCOL];
			char denc[MAXCOL];
			double tc = c(j);
			rat_approx(tc,numc,denc);
			atomp(numc,Num[j]);
			atomp(denc,Den[j]);
		}
		lrs_set_obj_mp(P,Q,Num,Den,MINIMIZE);
		lrs_clear_mp_vector(Num,d+1);
		lrs_clear_mp_vector(Den,d+1);
	};

	void makeH (lrs_dic *P, lrs_dat *Q,const mat<T>& M)
	{
		long row, j;
		long m=Q->m;       /* number of inequalities      */
		long n=Q->n;       /* hypercube has dimension n-1 */
		long d=P->d;
		lrs_mp_vector Num, Den;
		for (row=1;row<=m;row++){
			Num=lrs_alloc_mp_vector(d+1);
			Den=lrs_alloc_mp_vector(d+1);
			for(j=0;j<=d;j++){
				char num[MAXCOL];
				char den[MAXCOL];
				double t = M(row-1,j);
				rat_approx(t,num,den);
				atomp(num,Num[j]);
				atomp(den,Den[j]);
			}
			lrs_set_row_mp(P,Q,row,Num,Den,GE);
			lrs_clear_mp_vector(Num,d+1);
			lrs_clear_mp_vector(Den,d+1);
		}
	};

	/* generate H-representation of a unit hypercube */
	void makecube (lrs_dic *P, lrs_dat *Q)
	{
		long num[MAXCOL];
		long den[MAXCOL];
		long row, j;
		long m=Q->m;       /* number of inequalities      */
		long n=Q->n;       /* hypercube has dimension n-1 */
		for (row=1;row<=m;row++)
		{
			for(j=0;j<n;j++)
			{ num [j] = 0;
				den [j] = 1;
			}
			if (row < n)
			{ num[0] = 1;
				num[row] = -1;
			}
			else
			{ num[0] = 0;
				num[row+1-n] = 1;
			}
			lrs_set_row(P,Q,row,num,den,GE);
		}
		/* set up some objective function */
		 printf("\n\nObjective function: ");

		 for(j=0;j<n;j++)
		    { num [j] = 0;
		      den [j] = 1;
		      lprat(" ",num[j],den[j]);
		    }
		  lrs_set_obj(P,Q,num,den,MAXIMIZE);
	};

	/* print a row of A matrix in output in "original" form  */
	/* rowd+1 is the dimension of output vector                */
	/* if input is H-rep. output[0] contains the RHS      */
	/* if input is V-rep. vertices are scaled by 1/output[1] */
	// void lrs_print(mat<T>* o,long r, const char *name, lrs_dat * Q, lrs_mp_vector output, long rowd)
	// {
	//
	// 	for (long i = 0; i <= rowd; i++){
	// 		double d = mptoi(output[i]);
	// 		(*o)(r,i) = T(d);
	// 	}
	//
	// };
	mat<T>  ineq_print(lrs_mp_matrix Ain,lrs_dic *P,lrs_dat *Q)
	{
		long I;
		long *redineq=Q->redineq;
		long nlinearity = Q->nlinearity; /* number of linearities   */
		long m = P->m_A;
		/* restore as mplrs loses this */
		for (I = 0; I < nlinearity; I++){
			redineq[Q->linearity[I]]=2;
		}

		long nredund = 0;		/* count number of non-redundant inequalities */

		for (I = 1; I <= m; I++){
			if (redineq[I] == 0){
				nredund++;
			}
		}

		// fprintf (lrs_ofp, "\nbegin");
		// fprintf (lrs_ofp, "\n%ld %ld rational", nlinearity+nredund, Q->n);

		long C =0;
		mat<T> out(nredund,Q->inputd+1);
		for (I = 1; I <= m; I++){
			if (redineq[I] == 0){
				//lrs_print(&out,c,"", Q, Ain(i), Q->inputd);
				for (long X = 0; X < out.size_x; X++){
					out(C,X) = (*this)(I-1,X);
				}
				C++;
			}
		}

		return out;
	};
	vec<T> vert_print(lrs_dat * Q, lrs_mp_vector output)
	{
		vec<T> v(Q->n);
		for (long I = 1; I < Q->n; I++){
			double d;
			rattodouble(output[I],output[0],&d);
			v(I) = T(d);
		}
		return v;
	};
	mat<T> redund()
	{
		lrs_mp_matrix Ain;            /* holds a copy of the input matrix to output at the end */

		long ineq;			/* input inequality number of current index             */
		long *redineq;
		long i, j, d, m;
		long nlinearity;		/* number of linearities */
		long nredund;       /* number of redundant rows*/
		long lastdv;
		long index;			/* basic index for redundancy test */

		lrs_dic *P;	/* structure for holding current dictionary and indices  */
		lrs_dat *Q;	/* structure for holding static problem data             */
		lrs_mp_vector output;	/* one line of output:ray,vertex,facet,linearity */
		lrs_mp_matrix Lin;    /* holds input linearities if any are found      */

		/* Global initialization - done once */
		if ( !lrs_init ("")){
			cout << "lrs not initing" << endl;
			return mat<T>(1,1);
		}

		/* compute the vertices of a set of hypercubes given by */
		/* their H-representations.                             */
		/* allocate and init structure for static problem data */

		Q = lrs_alloc_dat ("");
		if (Q == NULL){
			cout << "lrs not allocing dat" << endl;
			return mat<T>(1,1);
		}

		/* now flags in lrs_dat can be set */

		Q->n=long(size_x);           /* number of input columns         (dimension + 1 )  */
		Q->m=long(size_y);         /* number of input rows = number of inequalities     */

		output = lrs_alloc_mp_vector (Q->n);

		P = lrs_alloc_dic (Q);   /* allocate and initialize lrs_dic      */
		if (P == NULL){
			cout << "lrs not allocing dic" << endl;
			return mat<T>(1,1);
		}

		/* Build polyhedron: constraints and objective */
		makeH(P,Q,(*this));

		/* if non-negative flag is set, non-negative constraints are not input */
		/* explicitly, and are not checked for redundancy                      */

		m = P->m_A;              /* number of rows of A matrix */
		d = P->d;
		redineq = Q->redineq;

		Q->Ain = lrs_alloc_mp_matrix (m, d);     /* make a copy of A matrix for output later            */
		Ain=Q->Ain;

		for (i = 1; i <= m; i++){
			for (j = 0; j <= d; j++){
				copy (Ain[i][j], P->A[i][j]);
			}
		}

		/*********************************************************************************/
		/* Step 1: Find a starting cobasis from default of specified order               */
		/*         Lin is created if necessary to hold linearity space                   */
		/*********************************************************************************/

		if (!lrs_getfirstbasis (&P, Q, &Lin, TRUE)){
			lrs_clear_mp_matrix(Ain,P->m_A,P->d);
			lrs_free_dic(P,Q);
			lrs_clear_mp_vector (output, Q->n);
			lrs_free_dat (Q);             /* deallocate lrs_dat */
			return mat<T>(1,1);
		}

		/* Pivot to a starting dictionary                      */
	  /* There may have been column redundancy               */
	  /* If so the linearity space is obtained and redundant */
	  /* columns are removed. User can access linearity space */
	  /* from lrs_mp_matrix Lin dimensions nredundcol x d+1  */

		/*********************************************************************************/
		/* Step 2: Test rows i where redineq[i]=1 for redundancy                         */
		/*********************************************************************************/

		/* note some of these may have been changed in getting initial dictionary        */

		m = P->m_A;
		d = P->d;
		nlinearity = Q->nlinearity;
		lastdv = Q->lastdv;

		/* linearities are not considered for redundancy */

		for (i = 0; i < nlinearity; i++){
			redineq[Q->linearity[i]] = 2L;
		}

		/* Q->verifyredund always false in lrs, set by mplrs to check duplicated redundancy removal */
		/* Q->noredundcheck overides this to skip verification                                      */

		if(Q->noredundcheck && Q->verifyredund){
			goto done;
		}

		/* mplrs sets redineq[i]==-1 for guaranteed redundant inequalities */
		/* these rows must be zeroed out before testing the others         */

		/* this is never run by lrs, final step of mplrs redund */
		if (Q->verifyredund){
			for (index = lastdv + Q->redineq[0]; index <= m + d; index++){
				ineq = Q->inequality[index - lastdv];     /* the input inequality number corr. to this index */
				if( redineq[ineq]== -1 ){
					checkindex (P, Q, -index);             /* used to zero correct row of A no LP solved */
				}
			}
		}

		/* rows 0..lastdv are cost, decision variables, or linearities  */
		/* other rows need to be tested                                */

		for (index = lastdv + Q->redineq[0]; index <= m + d; index++){
			ineq = Q->inequality[index - lastdv];	/* the input inequality number corr. to this index */
			Q->redineq[0] = ineq;                     /* used for restarting after arithmetic overflow    */
			if( redineq[ineq]==1 ){
				redineq[ineq] = checkindex (P, Q, index);
				if(!Q->mplrs && Q->verbose){
					if( redineq[ineq]==1 ){
						//lrs_printrow ("*re ", Q, Ain[ineq], Q->inputd);
					}
					else{
						//lrs_printrow ("*nr ", Q, Ain[ineq], Q->inputd);
					}
				}
			}
		}


		/* rows 0..lastdv are cost, decision variables, or linearities  */
		/* other rows need to be tested                                */

		done:

		m = P->m_A;              /* number of rows of A matrix */
		nlinearity = Q->nlinearity;
		/* restore as mplrs loses this */
		for (i = 0; i < nlinearity; i++){
			redineq[Q->linearity[i]]=2;
		}
		nredund = 0;		/* count number of non-redundant inequalities */
		for (i = 1; i <= m; i++){
			if (redineq[i] == 0){
				nredund++;
			}
		}

		mat<T> out = ineq_print(Ain,P,Q);

		lrs_clear_mp_matrix(Ain,P->m_A,P->d);

		//Q->m=P->m;

		lrs_free_dic(P,Q);

		lrs_clear_mp_vector (output, Q->n);

		lrs_free_dat (Q);             /* deallocate lrs_dat */
		//lrs_close ("");
		return out;
	};
	//tweaking some code from lrslib.h
	vec<T> lpsolve(const vec<T>& obj)
	{
		lrs_dic *P;	/* structure for holding current dictionary and indices  */
		lrs_dat *Q;	/* structure for holding static problem data             */
		lrs_mp_vector output;	/* one line of output:ray,vertex,facet,linearity */
		long col;	/* output column index for dictionary            */
		if ( !lrs_init ("")){return vec<T>(1);}
		Q = lrs_alloc_dat ("LRS globals");
		if (Q == NULL){return vec<T>(1);}
		/* now flags in lrs_dat can be set */
		//strcpy(Q->fname,"");
		Q->n=long(size_x);           /* number of input columns         (dimension + 1 )  */
		Q->m=long(size_y);         /* number of input rows = number of inequalities     */
		Q->lponly=TRUE;      /* we do not want all vertices generated!   */

		output = lrs_alloc_mp_vector (Q->n);

		P = lrs_alloc_dic (Q);   /* allocate and initialize lrs_dic      */

		if (P == NULL){return vec<T>(1);}

		/* Build polyhedron: constraints and objective */

		makeLP(P,Q,(*this), obj);

		//makecube(P,Q);

		/* Solve the LP  */
		if (!lrs_solve_lp(P,Q)){return vec<T>(1);}
		/* Print output */
		//prat ("\nObjective value = ", P->objnum, P->objden);
		double d;
		vec<T> v(Q->n);
		for (col = 0; col <= P->d; col++){
			if (lrs_getsolution(P, Q, output, col)){
				v = vert_print(Q, output);
				//lrs_printoutput (Q, output);
			}
		}
		rattodouble(P->objnum,P->objden,&d);
		v(0) = d;
		/* free space : do not change order of next lines! */
		lrs_clear_mp_vector (output, Q->n);
		lrs_free_dic (P,Q);         /* deallocate lrs_dic */
		lrs_free_dat (Q);           /* deallocate lrs_dat */
		//lrs_close("");
		return v;
	};
	//writes matrix to specified file
	void print(const std::string& Aname)
	{
		std::size_t X,Y;
		ofstream OSTRM(Aname);
		if (!OSTRM.is_open()) {
			cerr << "Output file isn't open" << endl;
			exit(1);
		}
		for ( Y = 0; Y < size_y; Y++){
			string sep = "";
			for (X = 0; X < size_x; X++){
				OSTRM << sep << A[Y+size_y*X];
				sep = ",";
			}
			OSTRM << endl;
		}
		OSTRM.close();
	}
	//writes matrix to .ine file for lrs
	void printIne(const std::string& Aname)
	{
		mat<T> temp(size_y,size_x);
		temp.set_all(*this);
		std::size_t X,Y;
		string sep = " ";
		ofstream OSTRM("temp.float");
		if (!OSTRM.is_open()) {
			cerr << "Output file isn't open" << endl;
			exit(1);
		}
		OSTRM <<  "frame" << endl;
		OSTRM <<  "*tempframe" << endl;
		OSTRM <<  "H-representation" << endl;
		OSTRM <<  "begin" << endl;
		OSTRM << size_y << " " << size_x << " rational"  << endl;
		for ( Y = 0; Y < size_y; Y++){
			for (X = 0; X < size_x; X++){
				OSTRM << temp(Y,X) << " ";
			}
			OSTRM << endl;
		}
		OSTRM << "end" << endl;
		OSTRM.close();
		f2r("temp.float",Aname.c_str());
	};
	//writes matrix as optimization .ine for lrs
	//maxmin: maximization or minimization bool
	//obj: the objective function
	void printLP(bool maxmin, const vec<T>& obj, const std::string& Aname)
	{
		mat<T> temp(*this);
		std::size_t X,Y;
		string sep = " ";
		ofstream OSTRM("templp.float");
		if (!OSTRM.is_open()) {
			cerr << "Output file isn't open" << endl;
			exit(1);
		}
		OSTRM <<  "frame" << endl;
		OSTRM <<  "*tempframe" << endl;
		OSTRM <<  "H-representation" << endl;
		OSTRM <<  "begin" << endl;
		OSTRM << size_y << " " << size_x << " rational"  << endl;
		for ( Y = 0; Y < size_y; Y++){
			for (X = 0; X < size_x; X++){
				OSTRM << temp(Y,X) << " ";
			}
			OSTRM << endl;
		}
		OSTRM << "end" << endl;
		OSTRM << "lponly" << endl;
		if (maxmin){
			OSTRM << "maximize" << endl;
		}else{
			OSTRM << "maximize" << endl;
		}
		for (X = 0; X < size_x; X++){
			OSTRM << obj(X) << " ";
		}
		OSTRM << endl;
		OSTRM.close();
		f2r("templp.float",Aname.c_str());
	};

	//writes to file with specified precision
	void print(const std::string& Aname,std::size_t precision)
	{
		std::size_t X,Y;
		ofstream OSTRM(Aname);
		OSTRM.precision(precision);
		OSTRM.setf(ios::fixed);
		if (!OSTRM.is_open()) {
			cerr << "Output file isn't open" << endl;
			exit(1);
		}
		for ( Y = 0; Y < size_y; Y++){
			string sep = "";
			for (X = 0; X < size_x; X++){
				OSTRM << sep << A[Y+size_y*X];
				sep = ",";
			}
			OSTRM << endl;
		}
		OSTRM.close();
	};
	//appends to specified file
	void printAppend(const std::string& Aname)
	{
		std::size_t X,Y;
		ofstream OSTRM;
		OSTRM.open(Aname, std::ios_base::app);
		if (!OSTRM.is_open()) {
			cerr << "Output file isn't open" << endl;
			exit(1);
		}
		for ( Y = 0; Y < size_y; Y++){
			string sep = "";
			for (X = 0; X < size_x; X++){
				OSTRM << sep << A[Y+size_y*X];
				sep = ",";
			}
			OSTRM << endl;
		}
		OSTRM.close();
	};
	//appends to specified file with select precision
	void printAppend(const std::string& Aname,std::size_t precision)
	{
		std::size_t X,Y;
		ofstream OSTRM;
		OSTRM.open(Aname, std::ios_base::app);
		if (!OSTRM.is_open()) {
			cerr << "Output file isn't open" << endl;
			exit(1);
		}
		for ( Y = 0; Y < size_y; Y++){
			string sep = "";
			for (X = 0; X < size_x; X++){
				OSTRM << sep << A[Y+size_y*X];
				sep = ",";
			}
			OSTRM << endl;
		}
		OSTRM.close();
	};
	//cout the matrix
	void print()
	{
		#ifdef DEBINFO
		cout << "mat[" << size_y<<"," << size_x << "]"<< endl;
		#endif
		for (std::size_t Y = 0; Y < size_y; Y++){
			for (std::size_t X = 0; X < size_x; X++){
				cout<<setprecision(3) << setfill(' ') << setw(9) << A[Y+size_y*X] << " ";
			}
			cout << endl;
		}
	};
};
//convert vector into a matrix 1 * size
template <typename T>
mat<T> to_mat(const vec<T>& v)
{
	mat<T> out(1,v.size);
	for (std::size_t X = 0; X < out.size_x; X++){
		out(0,X)=v(X);
	}
	return out;
};

template <typename T>
mat<T> make_diag(const vec<T>& v)
{
	mat<T> o(v.size,v.size);
	for ( std::size_t Y = 0; Y < v.size; Y++ ){
		o(Y,Y) = v(Y);
	}
	return o;
};

//reads lrs ine and returns a mat<T>
template <typename T>
mat<T> from_ine(const std::string& Iname)
{
	T i, j;
	std::size_t size_y = ineRow(Iname);
	std::size_t size_x = ineCol(Iname);
	size_t pos;
	mat<T> out(size_y,size_x);
	string line,str;
	T f;
	ifstream Ifile(Iname);
	istringstream issline, iss;
	std::size_t Y = 0;
	while(getline(Ifile,line) && !line.empty()) {
		if (line == "begin") break;
	}
	getline(Ifile,line);
	while(getline(Ifile,line) && !line.empty()) {
		if (line == "end") break;
		istringstream issline(line);
		std::size_t X = 0;
		while(getline(issline >> std::ws,str,' ')){
			pos = str.find('/');// fraction
			if (pos != string::npos) {
				// get numerator
				iss.str(str.substr(0, pos));
				iss >> i;
				iss.clear();
				// get denominator
				iss.str(str.substr(pos + 1));
				iss >> j;
				iss.clear();
				// calculate fraction
				f = (double) i / (double) j;
			}
			else{
				// otherwise convert int to double
				f = strtof(str.c_str(), NULL);
			}
			out(Y,X) = f;
			X++;
			str.clear();
		}
		issline.clear();
		Y++;
	}
	Ifile.close();
	return out;
};

#endif
