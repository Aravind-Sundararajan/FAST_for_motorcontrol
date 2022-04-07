#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
int main(int argc, char *argv[])
{
  if (argc != 4) {
    cout << "Usage: ./redunder in out is_ine_bool" << endl;
    exit(0);
  }
  int type = atoi(argv[3]);
  if (type == 0){
    mat<double> M(argv[1]);
    mat<double> O = M.redund();
    O.print(argv[2]);
  }
  else{
    mat<double> M = from_ine<double>(argv[1]);
    mat<double> O = M.redund();
    O.print(argv[2]);
  }
  return 0;
}
