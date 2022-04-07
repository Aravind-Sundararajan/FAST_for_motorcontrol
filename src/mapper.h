#ifndef MAPPER_H_
#define MAPPER_H_
#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"

using namespace std;
int map(const std::string& Aname,const std::string& bname,const std::string& out_name)
{
  mat<double> A(Aname);
  ofstream OSTRM(out_name);
  if (!OSTRM.is_open()) {
    cerr << "Output file isn't open" << endl;
    exit(1);
  }
  string str,line;
  double f;
  size_t sz;
  ifstream bfile(bname);
  std::size_t Y = 0;
  vec<double> b(A.size_x);
  //cout << endl;
  while(getline(bfile,line) && !line.empty()) {
    istringstream issline(line);
    std::size_t X = 0;
    while(getline(issline,str,',')){
      f = stof(str,&sz);
      if (abs(f) < TOL){
        f=0;
      }
      if (X>0){
        b(X-1) = f;
      }
      X++;
    }
    b(X-1) =1;
    issline.clear();
    string sep = "";
    for (std::size_t Y = 0; Y < A.size_y; Y++){
      double vec=0;
      for (std::size_t X = 0; X < A.size_x; X++){
        vec+=b(X)*A(Y,X);
      }
      OSTRM << sep << vec;
      sep=",";
    }
    Y++;
    OSTRM << endl;
  }
  OSTRM.close();
  return 0;
};
#endif
