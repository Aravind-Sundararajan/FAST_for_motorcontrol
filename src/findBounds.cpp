#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
using namespace std;
typedef std::numeric_limits< double > dbl;
std::size_t nCol(char* fname)
{
std::size_t nCols = 0;
string line,temp;
ifstream file(fname);
getline(file,line);
istringstream issline(line);
while(getline(issline,temp,',')){
nCols++;
}
return nCols;
}
int main(int argc, char **argv)
{
string line, str, sep;
double f;
  if (argc != 3) {
    cout << "Usage: ./findBounds input-file output-file" << endl;
    exit(0);
  }
#ifdef TIME_PROFILE
  auto start = chrono::high_resolution_clock::now();
#endif
  size_t sz;
  std::size_t nc = nCol(argv[1]);
  cout << "there were: " << to_string(nc) << " columns." << endl;
  mat<double> MinMax(2,nc);
  MinMax.set_all((double)0);
  ifstream file(argv[1]);
  std::size_t rid = 0;
  while(getline(file,line) && !line.empty()) {
    istringstream issline(line);
    std::size_t cid = 0;
    while(getline(issline,str,',')){
      f = stof(str,&sz);
      if (rid == 0){
        MinMax(0,cid) = f;//max
        MinMax(1,cid) =  f;//min
      }
      if (MinMax(0,cid)<f){
          MinMax(0,cid)=f;
      }
      if (MinMax(1,cid)>f){
	      MinMax(1,cid)=f;
      }
	cid++;
    }
	rid++;
  }
  MinMax.print(argv[2]);
#ifdef TIME_PROFILE
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << endl << "Finding the bounds took "
       << duration.count() << " microseconds." << endl;
#endif
  return 0;
}
