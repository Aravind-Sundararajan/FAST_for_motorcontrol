#ifndef UTILS_H_
#define UTILS_H_
#include "base.h"
using namespace std;
//returns number of columns from input matrix file
inline std::size_t nCol( const std::string& Aname)
{
	string str,line;
	ifstream Afile(Aname);
	std::size_t X = 0;
	getline(Afile,line);
	istringstream issline(line);
	while(getline(issline,str,',')){
		X++;
	}
	issline.clear();
	Afile.close();
	return X;
};
//returns number of rows in input matrix file
inline std::size_t nRow( const std::string& Aname)
{
	string str,line;
	ifstream Afile(Aname);
	std::size_t Y = 0;
	while(getline(Afile,line) && !line.empty()) {
		istringstream issline(line);
		Y++;
		issline.clear();
	}
	Afile.close();
	return Y;
};
//return number of rows from an input lrs .ine file
inline std::size_t ineRow(const std::string& Iname)
{
	string line,str;
	ifstream Ifile(Iname);
	while(getline(Ifile,line) && !line.empty()) {
		if (line == "begin") break;
	}
	std::size_t Y = 0;
	getline(Ifile,line);
	while(getline(Ifile,line) && !line.empty()) {
		if (line == "end") break;
		Y++;
	}
	Ifile.close();
	return Y;
};
//return number of cols from an input lrs .ine file
inline std::size_t ineCol(const std::string& Iname)
{
	string line,str;
	ifstream Ifile(Iname);
	while(getline(Ifile,line) && !line.empty()) {
		if (line == "begin") break;
	}
	getline(Ifile,line);
	getline(Ifile,line);
	istringstream issline(line);
	std::size_t X = 0;
	while(getline(issline >> std::ws,str,' ')){
		X++;
	}
	Ifile.close();
	return X;
};
//return string rational approximation
inline int f2r(const char* Iname,const char* Oname)
{
	long int  m,n;
	int i,j;
	long atol();
	FILE *fin;/* input file pointer       */
	FILE *fout;/* output file pointer       */
	char  buf[BUFSIZ];
	fin = fopen (Iname, "r+");
	fout = fopen (Oname, "w+");
	while ( fgets(buf,BUFSIZ,fin) !=NULL ){
		fputs(buf,fout);
		if (strncmp(buf,"begin",5)==0) break;
	}
	if (fscanf(fin,"%ld %ld %s",&m,&n,buf)==EOF){
		fprintf(stderr,"No begin line");
		return 1;
	}
	fprintf(fout,"%ld %ld rational\n",m,n);
	for (i=0;i<m;i++)   {
		for(j=0;j<n;j++)	{
			char *p;
			char *frac;
			std::size_t k;
			fscanf(fin,"%s",buf);
			if ((p=strchr(buf,'.'))){
				*p=0;
				frac=&p[1];
				fprintf(fout,"%s%s/1",buf,frac);
				for (k=0; k<strlen(frac); k++)
				fprintf(fout,"0");
			}else {
				fprintf(fout,"%s",buf);
			}
			fprintf(fout," ");
		}
		fputs("\n",fout);
	}
	fgets(buf,BUFSIZ,fin);
	while (fgets(buf,BUFSIZ,fin) !=NULL ) {
		fputs(buf,fout);
		if (strncmp(buf,"maximize",8)==0) break;
		if (strncmp(buf,"minimize",8)==0) break;
	}
	for(j=0;j<n;j++)	{
		char *p;
		char *frac;
		std::size_t k;
		fscanf(fin,"%s",buf);
		if ((p=strchr(buf,'.'))){
			*p=0;
			frac=&p[1];
			fprintf(fout,"%s%s/1",buf,frac);
			for (k=0; k<strlen(frac); k++)
			fprintf(fout,"0");
		}else {
			fprintf(fout,"%s",buf);
		}
		fprintf(fout," ");
	}
	fputs("\n",fout);
	fclose(fin);
	fclose(fout);
	free(fin);
	free(fout);
	return 0;
};
//return string rational approximation
void rat_approx(double v,char num[], char den[])
{
	char  buf[BUFSIZ];
	sprintf(buf,"%.17f",v);
	char *p;
	char *frac;
	std::size_t k;
	if ((p=strchr(buf,'.'))){
		*p=0;
		frac=&p[1];
		sprintf(num,"%s%s",buf,frac);
		sprintf(den,"1");
		for (k=0; k<strlen(frac); k++)
		sprintf(den+strlen(den),"%s","0");
	}
	else {
		sprintf(num,"%s",buf);
		sprintf(den,"%s","1");
	}
};
//some convenience functions for mat.fun(f) or vec.fun(f)
template <typename T>
T powa(T n,double d){return pow(n,d);}
template <typename T>
T pow2(T n){return powa(n,2);}
template <typename T>
T powroot(T n){return powa(n,.5);}
template <typename T>
T recip(T n){return 1/n;};
//returns a string "true" or "false" instead of bool 0 1
template <typename T>
const char* boolstring(bool b)
{
	return b ? "true" : "false";
};
#endif
