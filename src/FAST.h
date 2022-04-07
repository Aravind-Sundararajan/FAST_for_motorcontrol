#ifndef FAST_H_
#define FAST_H_
#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
#include "polytope.h"
#include "optim.h"
#include "redund.h"
#include "window.h"
#include <sys/stat.h>
using namespace std;
int run_fast(int F, double DT, double TCA, double TCD, const std::string& path, const std::string& outdir, int niter)
{
  std::ostringstream ss;
  for (std::size_t n = 0; n<niter; n++){
    cout <<"deltas Path:"<< path+"deltas.csv" << endl;
    mat<double> skew(path+"deltas.csv");
    window<double> w(F,DT,TCA,TCD,path);
    w.run();
    ss.str("");
    ss << std::setfill('0') << std::setw(3) << n;
    w.printtraj(outdir+ss.str() + ".csv");
    // w.printlb(outdir+ss.str() + ".lb");
    // w.printub(outdir+ss.str() + ".ub");
  }
  return 0;
}
#endif
