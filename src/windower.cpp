#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
#include "polytope.h"
#include "optim.h"
#include "redund.h"
#include "window.h"
#include "FAST.h"
int main(int argc, char *argv[])
{//fingerprints of the gods
  if (argc != 8) {
    cout << "Usage: ./windower nframes DT TCA TCD path_in path_out niter" << endl;
    exit(0);
  }
  std::size_t F = atoi(argv[1]);
  double DT = atof(argv[2]);
  double TCA = atof(argv[3]);
  double TCD = atof(argv[4]);
  std::string indir = argv[5];
  std::string outdir = argv[6];
  std::size_t niter = atoi(argv[7]);
  std::ostringstream ss;
  run_fast(F,DT,TCA, TCD, indir, outdir, niter);
}
