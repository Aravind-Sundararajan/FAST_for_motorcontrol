#include "mapper.h"
using namespace std;
int main(int argc, char **argv)
{
	cout << argc << endl;
	if (argc != 6) {
		cout << "Usage: ./fsm mapfile spacefile outfile" << endl;
		exit(0);
	}
	map(argv[1],argv[2],argv[3]);
}
