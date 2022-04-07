#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
#include "polytope.h"
#include "optim.h"
#include "redund.h"
#include "window.h"
#include <sys/stat.h>

int main(int argc, char *argv[])
{
  if (argc != 7) {
    cout << "Usage: ./windower nframes DT TCA TCD path_in path_out" << endl;
    exit(0);
  }
  std::size_t F = atoi(argv[1]);
  double DT = atof(argv[2]);
  double TCA = atof(argv[3]);
  double TCD = atof(argv[4]);
  std::string path = argv[5];
  std::string outdir = argv[6];
  mat<double> H(path+"000.H");
  mat<double> skew(path+"delta.csv");
  std::size_t n_muscles = H.size_x;
  std::ostringstream ss;

  window<double> wub(F,DT,TCA,TCD,path);
  wub.run_optim_all(-1.0);
  wub.printtraj(outdir + "maxall.csv");

  window<double> wlb(F,DT,TCA,TCD,path);
  wlb.run_optim_all(1.0);
  wlb.printtraj(outdir + "minall.csv");

  for (std::size_t m = 0; m<n_muscles; m++){
    ss.str("");
    ss << std::setfill('0') << std::setw(3) << m;

    window<double> w(F,DT,TCA,TCD,path);
    w.run_optim_one(m,1.0);
    w.printtraj(outdir+ss.str() + "max.csv");

    window<double> w2(F,DT,TCA,TCD,path);
    w2.run_optim_one(m,-1.0);
    w2.printtraj(outdir+ss.str() + "min.csv");

  }



  return 0;
}
