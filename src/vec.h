#ifndef VEC_H_
#define VEC_H_
#include "base.h"
#include "util.h"
//container template class that holds a vector as a 1D array
//convenient operators, arbitrary function application on elem
//rowops, printers, etc...
using namespace std;
template <typename T>
class vec
{
public:
	T *A{};
	std::size_t size{};
	//default constructor
	vec()
	{
		size = (std::size_t)1;
		A = new T [size]{};
		*(A) = (T)0;
	};
	//size constructor
	vec(std::size_t s)
	{
		size = s;
		A = new T [size]{};
	};
	//finput constructor
	vec(const std::string& Aname)
	{
		size = nElem(Aname);
		A = new T [size]{};
		set_all(Aname);
	};
	//copy constructor
	vec(const vec<T>& v)
	{
		size = v.size;
		A = new T [size]{};
		for (std::size_t X =0; X< size; X++){
			A[X] = v(X);
		}
	};
	//copy assignment
	vec<T>& operator=(const vec<T>& v){
		size = v.size;
		delete[] A;
		A = new T [size]{};
		set_all(v);
		return *this;
	}
	//destructor
	~vec()
	{
		delete[] A;
	};
	//idk why I have this method
	void Erase()
	{
		delete[] A;
		A = nullptr;
		size = 0;
	};
	//template cast assignment
	template <typename U>
	vec<T>& operator=(const vec<U>& v)
	{
		delete[] A;
		size = v.size;
		A = new T [size]{};
		set_all(v);
		return *this;
	};
	//template cast constructor
	template <typename U>
	vec(const vec<U>& v)
	{
		size = v.size;
		A = new T [size]{};
		for (std::size_t X =0; X< size; X++){
			A[X] = v(X);
		}
	}
	//return value at X
	T& operator()(std::size_t X)
	{
		return A[X];
	};
	//return const value at X
	T& operator()(std::size_t X) const
	{
		return A[X];
	};
	//funny looking assignment, kinda like a tuple
	void operator()(std::size_t X, T v)
	{
		A[X] = v;
	};
	//return value at X
	T& get(std::size_t X)
	{
		return A[X];
	};
	//return const value at X
	T& get(std::size_t X) const
	{
		return A[X];
	};
	//assign value v to position X
	void set(std::size_t X,T v)
	{
		A[X] = v;
	};
	//set every index to v
	void set_all(T v)
	{
		for (std::size_t X = 0; X < size; X++){
			set(X, v);
		}
	};
	//copy assign each index
	void set_all(const vec<T>& v)
	{
		std::size_t sz;
		if (v.size < size){
			sz = v.size;
		}else {
			sz = size;
		}
		for (std::size_t X = 0; X < sz; X++){
			set(X, v(X));
		}
	}
	// set values from file
	void set_all(const std::string& Aname)
	{
		string line,str;
		T f;
		ifstream Afile(Aname);
		size_t X = 0;
		while(getline(Afile,line) && !line.empty()) {
			istringstream issline(line);
			while(getline(issline,str,',')){
				istringstream ss(str);
				ss >> f;
				A[X]=f;
				X++;
				ss.clear();
			}
			issline.clear();
		}
		Afile.close();
	};
	//append scalar value to end of vector
	vec<T> append (T v)
	{
		vec<T> u(size+1);
		u.set_all(*this);
		u(size) = v;
		return u;
	};
	//remove value at index from vector and return new vectors size-1
	vec<T> cutrow(std::size_t Y)
	{
		vec<T> out(1);
		for (std::size_t YY =0; YY< size; YY++){
			if (YY!=Y){
				out=out.append((*this)(YY));
			}
		}
		out=out.trunc(1,size-1);
		return out;
	};
	//appends vector to end of vector
	vec<T> concat(const vec<T>& b)  const
	{
		std::size_t X;
		vec<T> out(size);
		out.set_all(*this);
		for (X =0; X <b.size; X++){
			out = out.append(b(X));
		}
		return out;
	};
	//returns sum of all values in vector
	T sum()
	{
		T out = 0;
		std::size_t X;
		for (X =0; X< size; X++){
			out = out + get(X);
		}
		return out;
	};
	//arbitrary function application on each elem of vector
	vec<T> fun(T f(T)) // pass function pointer
	{
		vec<T> out(size);
		std::size_t X;
		for (X =0; X< size; X++){
			out(X)=f((*this)(X));
		}
		return out;
	};
	//truncate vector starting at position X0 to size szx
	vec<T> trunc (std::size_t X0,std::size_t szx)
	{
		vec<T> out(szx);
		for (std::size_t X =X0; X< X0+szx; X++){
			out(X-X0,get(X));
		}
		return out;
	};
	//OPERATORS
	vec<T> operator+ ()
	{
		vec<T> u(*this);
		return u;
	};
	vec<T> operator- ()
	{
		vec<T> u(*this);
		u*=(T)-1;
		return u;
	};
	vec<T> operator+ (T v)
	{
		vec<T> u(*this);
		u+=v;
		return u;
	};
	vec<T> operator- (T v)
	{
		vec<T> u(*this);
		u-=v;
		return u;
	};
	vec<T> operator+ (const vec<T>& v) const
	{
		vec<T> u(*this);
		u+=v;
		return u;
	};
	vec<T> operator- (const vec<T>& v) const
	{
		vec<T> u(*this);
		u-=v;
		return u;
	};
	vec<T> operator* (T v)
	{
		vec<T> u(*this);
		u*=v;
		return u;
	};
	vec<T> operator/ (T v)
	{
		vec<T> u(*this);
		u/=v;
		return u;
	};
	vec<T> operator* (T v) const
	{
		vec<T> u(*this);
		u*=v;
		return u;
	};
	vec<T> operator/ (T v) const
	{
		vec<T> u(*this);
		u/=v;
		return u;
	};
	vec<T> operator* (const vec<T>& v)
	{
		vec<T> u(*this);
		u*=v;
		return u;
	};
	vec<T> operator/ (const vec<T>& v)
	{
		vec<T> u(*this);
		u/=v;
		return u;
	};
	T operator% (const vec<T>& v) const
	{
		T out = 0;
		for (std::size_t X = 0; X < size; X++){
			out+= A[X] * v(X);
		}
		return out;
	}
	//ASSIGNMENT OPERATORS
	vec<T>& operator+=(T v)
	{
		for (std::size_t X = 0; X < size; X++){
			A[X] = A[X] + v;
		}
		return *this;
	};
	vec<T>& operator-=(T v)
	{
		for (std::size_t X = 0; X < size; X++){
			A[X] = A[X] - v;
		}
		return *this;
	};
	vec<T>& operator+=(const vec<T>& v)
	{
		for (std::size_t X = 0; X < size; X++){
			A[X] = A[X] + v(X);
		}
		return *this;
	};
	vec<T>& operator-=(const vec<T>& v)
	{
		for (std::size_t X = 0; X < size; X++){
			A[X] = A[X] - v(X);
		}
		return *this;
	};
	vec<T>& operator*=(T v)
	{
		for (std::size_t X = 0; X < size; X++){
			A[X] = A[X] * v;
		}
		return *this;
	};
	vec<T>& operator/=(T v)
	{
		for (std::size_t X = 0; X < size; X++){
			A[X] = A[X] / v;
		}
		return *this;
	};
	vec<T>& operator*=(const vec<T>& v)
	{
		for (std::size_t X = 0; X < size; X++){
			A[X] = A[X] * v(X);
		}
		return *this;
	};
	vec<T>& operator/=(const vec<T>& v)
	{
		for (std::size_t X = 0; X < size; X++){
			A[X] = A[X] / v(X);
		}
		return *this;
	};
	//PRINTERS
	void print();
	void print(const std::string& Aname);
	void printAppend(const std::string& Aname);
	void print(const std::string& Aname,const std::size_t precision);
	void printAppend(const std::string& Aname,const std::size_t precision);
	std::size_t nElem(const std::string& Aname);
};
//print to cout
template <typename T>
void vec<T>::print()
{
	for (std::size_t X = 0; X < size; X++){
		cout <<setprecision(3) << setfill(' ') << setw(7) << get(X) << " ";
	}
	cout << endl;
};
//prints file
template <typename T>
void vec<T>::print(const std::string& Aname)
{
	std::size_t X;
	ofstream OSTRM(Aname);
	if (!OSTRM.is_open()) {
		cerr << "Output file isn't open" << endl;
		exit(1);
	}
	string sep = "";
	for (X = 0; X < size; X++){
		OSTRM << sep << get(X);
		sep = ",";
	}
	OSTRM << endl;
	OSTRM.close();
};
//prints to file at precision
template <typename T>
void vec<T>::print(const std::string& Aname,const std::size_t precision)
{
	std::size_t X;
	ofstream OSTRM(Aname);
	OSTRM.precision(dbl::max_digits10);
	OSTRM.setf(ios::fixed);
	if (!OSTRM.is_open()) {
		cerr << "Output file isn't open" << endl;
		exit(1);
	}
	string sep = "";
	for (X = 0; X < size; X++){
		OSTRM << sep << get(X);
		sep = ",";
	}
	OSTRM << endl;
	OSTRM.close();
};
//appends to end of file
template <typename T>
void vec<T>::printAppend(const std::string& Aname)
{
	std::size_t X;
	ofstream OSTRM;
	OSTRM.open(Aname, std::ios_base::app);
	if (!OSTRM.is_open()) {
		cerr << "Output file isn't open" << endl;
		exit(1);
	}
	string sep = "";
	for (X = 0; X < size; X++){
		OSTRM << sep << get(X);
		sep = ",";
	}
	OSTRM << endl;
	OSTRM.close();
};
//prints to end of file at precision
template <typename T>
void vec<T>::printAppend(const std::string& Aname,const std::size_t precision)
{
	std::size_t X;
	ofstream OSTRM;
	OSTRM.open(Aname, std::ios_base::app);
	OSTRM.precision(dbl::max_digits10);
	OSTRM.setf(ios::fixed);
	if (!OSTRM.is_open()) {
		cerr << "Output file isn't open" << endl;
		exit(1);
	}
	string sep = "";
	for (X = 0; X < size; X++){
		OSTRM << sep << get(X);
		sep = ",";
	}
	OSTRM << endl;
	OSTRM.close();
};
//counts number of elements in vec from input file
template <typename T>
std::size_t vec<T>::nElem(const std::string& Aname)
{
	string str,line;
	ifstream Afile(Aname);
	std::size_t X = 0;
	while(getline(Afile,line) && !line.empty()) {
		istringstream issline(line);
		while(getline(issline,str,',')){
			X++;
		}
		issline.clear();
	}
	Afile.close();
	return X;
};
//reads lrs lp .ext and returns vec
template <typename T>
inline vec<T> from_lp(const std::string& Iname)
{
	T i, j, f;
	size_t pos,X;
	vec<T> out(1);
	string line,str;
	ifstream Ifile(Iname);
	istringstream iss;
	while(getline(Ifile,line) && !line.empty()){
		if (line[0] == ' ') break;
	}
	istringstream issline(line);
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
		out = out.append(f);
		X++;
		str.clear();
	}
	Ifile.close();
	out = out.trunc(2,out.size-2);
	return out;
};
#endif
