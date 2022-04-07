#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
using namespace std;
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
		cout << "Usage: ./findMean input-file output-file" << endl;
		exit(0);
	}
	#ifdef TIME_PROFILE
	auto start = chrono::high_resolution_clock::now();
	#endif
	size_t sz;
	std::size_t nLines = 0;
	std::size_t nc = nCol(argv[1]);
	cout << "there were: " << to_string(nc) << " columns." << endl;
	vec<double> sums(nc);
	for (std::size_t id = 0; id < nc; id ++){
		sums(id) = 0;
	}
	ifstream file(argv[1]);
	while(getline(file,line) && !line.empty()) {
		nLines++;
		istringstream issline(line);
		std::size_t cid = 0;
		while(getline(issline,str,',')){
			f = stof(str,&sz);
			sums(cid)= sums(cid)+f;
			cid++;
		}
	}
	double nl = nLines;
	vec<double> Mean = sums/nl;
	Mean.print(argv[2]);
	cout << endl << "There were " << to_string(nLines) << " rows." << endl;
	#ifdef TIME_PROFILE
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	cout << endl << "Finding the mean took "
	<< duration.count() << " microseconds." << endl;
	#endif
	return 0;
}
