#include "base.h"
#include "util.h"
#include "kvp.h"
#include "vec.h"
#include "mat.h"
#include "polytope.h"
#include "hitandrun.h"

using namespace std;
int main(int argc, char **argv)
{
	cout << argc << endl;
	if (argc != 9) {
		cout << "Usage: ./sampler Afile bfile mfile outname n_samples thin" << endl;
		exit(0);
	}
	int n_samples = atoi(argv[4]);
	int thin      = atoi(argv[5]);
	polytope<double> P(argv[1],argv[2]);
	// P->print_A();
	// P->print_b();
	// P->print_H();
	vec<double> point(argv[3]);
	if (P.check_inside(point)){
		hitandrun<double> HAR(P, point, n_samples, thin);
	HAR.run(argv[4]);
	}else{
	cout << "cant run hit and run. The starting position wasn't inside H-rep." << endl;
	}
	return 1;
}
