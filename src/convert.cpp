#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
typedef std::numeric_limits< double > dbl;
using namespace std;
#ifdef USE_OUTPUT_FILE
#define OSTRM outfile
#else
#define OSTRM cout
#endif
int main(int argc, char **argv)
{
  istringstream issline, iss;
  string line, str;
  size_t pos;
  double i, j;
  double f;
  if (argc != 3) {
    cout << "Usage: ./convert input-file output-file" << endl;
    exit(0);
  }
  ifstream file(argv[1]);
#ifdef USE_OUTPUT_FILE
  ofstream outfile(argv[2]);
  if (!outfile.is_open()) {
    cerr << "Output file isn't open" << endl;
    exit(1);
  }
#endif
  OSTRM.precision(dbl::max_digits10);
  OSTRM.setf(ios::fixed);
#ifdef TIME_PROFILE
  auto start = chrono::high_resolution_clock::now();
#endif
  while (getline(file, line) && !line.empty()) {
    // ignore header and footer
    if (line[0] != ' ') continue;
    // convert line str to sstream
    istringstream issline(line);
	   string sep = "";
    // get elements from line
    while (issline >> str) {
      pos = str.find('/');
      if (pos != string::npos) {
        // fraction
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
        OSTRM.precision(dbl::max_digits10);
      } else {
        // otherwise convert int to double
        f = strtof(str.c_str(), NULL);
        OSTRM.precision(dbl::max_digits10);
      }
      OSTRM << sep << f;
	  sep = ",";
    }
    issline.clear();
    OSTRM << endl;
  }
OSTRM << endl;
#ifdef TIME_PROFILE
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << endl << "Converting from MPLRS rational format to decimal CSV took "
       << duration.count() << " microseconds." << endl;
#endif
#ifdef USE_OUTPUT_FILE
  outfile.close();
#endif
  return 0;
}
