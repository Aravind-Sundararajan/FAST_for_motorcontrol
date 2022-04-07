#include "base.h"
#include <boost/math/distributions.hpp>
using namespace std;
using namespace boost::math;
int nCols(char* input_name)
{
  int nCol = 0;
  string line,temp;
  ifstream file(input_name);
  getline(file,line);
  istringstream issline(line);
  while(getline(issline,temp,',')){
    nCol++;
  }
  issline.clear();
  file.close();
  return nCol;
}
int nRows(char* input_name)
{
  int nRow = 0;
  string line;
  ifstream file(input_name);
  while(getline(file,line) && !line.empty()){
    nRow++;
  }
  file.close();
  return nRow;
}
double* compute_mean(char* input_name, double mean[])
{
  int nCol = nCols(input_name);
  int nLines=0;
  double f;
  string line, str;
  size_t sz;
  //cout << "there were: " << to_string(nCol) << " columns." << endl;
  for (int id = 0; id < nCol; id ++){
  mean[id] = 0;
  }
  ifstream file(input_name);
  while(getline(file,line) && !line.empty()) {
  nLines++;
    istringstream issline(line);
    int cid = 0;
    while(getline(issline,str,',')){
      f = stof(str,&sz);
      mean[cid]+=f;
      cid++;
    }
    issline.clear();
  }
  file.close();
  for (int X = 0; X < nCol; X ++){
    mean[X] = mean[X]/nLines;
    //cout << mean[X] << " ";
  }
  //cout << endl;
  return mean;
}
double* compute_standard_deviation(char* input_name, double sd[])
{
  size_t sz;
  string line, str;
  int nCol = nCols(input_name);
  int nRow = nRows(input_name);
  double sos[nCol],m[nCol],f;
  for (int x = 0; x < nCol; x++){
    sos[x] = 0;
    m[x] = 0;
  }
  compute_mean(input_name,m);
  ifstream file(input_name);
  while(getline(file,line) && !line.empty()) {
    istringstream issline(line);
    int cid = 0;
    while(getline(issline,str,',')){
      f = stof(str,&sz);
      sos[cid]+= powf(f-m[cid],2);
    //  cout << sos[cid] << " ";
      cid++;
    }
    issline.clear();
  //  cout << endl;
  }
  file.close();
  for (int X = 0; X < nCol; X++){
    sd[X] = sqrt(sos[X]/(nRow-1));
  //  cout << sd[X] << " ";
  }
  //cout << endl;
  //cout << "nRow was: " << nRow <<". nCol was " << nCol << "." <<endl;
  return sd;
}
double* compute_CI(char* input_name,double CI[], double alpha)
{
  int nCol = nCols(input_name);
  int nRow = nRows(input_name);
  double mean[nCol], sd[nCol];
  compute_mean(input_name,mean);
  compute_standard_deviation(input_name,sd);
  normal_distribution<> dist;
  double Z = quantile(dist, 1-alpha / 2);
  //cout << "Z score is:" << T << endl;
  double w = 0;
  for (int X =0; X < nCol; X++){
    w = Z * sd[X] / sqrt(double(nRow));
    CI[0 + 2*X] = mean[X] - w;
    CI[1 + 2*X] = mean[X] + w;
    //cout << "{"<< CI[0 + 2*X] << "," <<  CI[1 + 2*X] << "}" << endl;
  }
  return CI;
}
void print_mean(char* mean_file,double mean[],int nCol)
{
ofstream OSTRM(mean_file);
OSTRM.precision(double_DECIMALS);
OSTRM.setf(ios::fixed);
string sep = "";
for (int X = 0;X < nCol; X++) {
  OSTRM << sep << mean[X];
  sep=",";
}
OSTRM << endl;
OSTRM.close();
}
void print_standard_deviation(char* sd_file,double sd[],int nCol)
{
  ofstream OSTRM(sd_file);
  OSTRM.precision(double_DECIMALS);
  OSTRM.setf(ios::fixed);
  string sep = "";
  for (int X = 0;X < nCol; X++) {
    OSTRM << sep << sd[X];
    sep=",";
  }
  OSTRM << endl;
  OSTRM.close();
}
void print_CI(char* CI_file,double CI[], int nCol)
{
  ofstream OSTRM(CI_file);
  OSTRM.precision(double_DECIMALS);
  OSTRM.setf(ios::fixed);
  for (int Y = 0; Y<2; Y++){
    string sep = "";
    for (int X = 0;X < nCol; X++){
      OSTRM << sep << CI[Y+ 2*X];
      sep=",";
    }
    OSTRM << endl;
  }
  OSTRM.close();
}
int main(int argc, char **argv)
{
if (argc != 6) {
  cout << "Usage: ./stats input-file output-mean-file output-sd-file output-CI-file alpha" << endl; //output-mean-file output-sd-file output-CI95-file
  exit(0);
}
int nCol = nCols(argv[1]);
//int nRow = nRows(argv[1]);
double alpha = (double)atoi(argv[5])/100;
//cout << "alpha is: " << alpha << endl;
double mean[nCol], sd[nCol], CI[2 * nCol];
compute_mean(argv[1],mean);
compute_standard_deviation(argv[1],sd);
compute_CI(argv[1],CI,alpha);
print_mean(argv[2],mean,nCol);
print_standard_deviation(argv[3],sd,nCol);
print_CI(argv[4],CI,nCol);
/*
#ifdef USE_OUTPUT_FILE
  ofstream outfile(argv[2]);
  if (!outfile.is_open()) {
    cerr << "Output file isn't open" << endl;
    exit(1);
  }
#endif
#ifdef TIME_PROFILE
  auto start = chrono::high_resolution_clock::now();
#endif
  // OSTRM.precision(double_DECIMALS);
  // OSTRM.setf(ios::fixed);
  // sep = "";
  // for (int id = 0;id < nc; id++) {
  //   OSTRM << sep << Mean[id];
  //   sep=",";
  // }
  // OSTRM << endl;
cout << endl << "There were " << to_string(nLines) << " rows." << endl;
#ifdef TIME_PROFILE
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << endl << "Finding the mean took "
       << duration.count() << " microseconds." << endl;
#endif
*/
return 0;
}
#include "base.h"
#include "util.h"
#include "kvp.h"
#include "vec.h"
#include "mat.h"
#include "polytope.h"
#include <boost/math/distributions.hpp>
using namespace std;
using namespace boost::math;
int nCols(char* input_name)
{
  int nCol = 0;
  string line,temp;
  ifstream file(input_name);
  getline(file,line);
  istringstream issline(line);
  while(getline(issline,temp,',')){
    nCol++;
  }
  issline.clear();
  file.close();
  return nCol;
}
int nRows(char* input_name)
{
  int nRow = 0;
  string line;
  ifstream file(input_name);
  while(getline(file,line) && !line.empty()){
    nRow++;
  }
  file.close();
  return nRow;
}
double* compute_mean(char* input_name, double mean[])
{
  int nCol = nCols(input_name);
  int nLines=0;
  double f;
  string line, str;
  size_t sz;
  //cout << "there were: " << to_string(nCol) << " columns." << endl;
  for (int id = 0; id < nCol; id ++){
  mean[id] = 0;
  }
  ifstream file(input_name);
  while(getline(file,line) && !line.empty()) {
  nLines++;
    istringstream issline(line);
    int cid = 0;
    while(getline(issline,str,',')){
      f = stof(str,&sz);
      mean[cid]+=f;
      cid++;
    }
    issline.clear();
  }
  file.close();
  for (int X = 0; X < nCol; X ++){
    mean[X] = mean[X]/nLines;
    //cout << mean[X] << " ";
  }
  //cout << endl;
  return mean;
}
double* compute_standard_deviation(char* input_name, double sd[])
{
  size_t sz;
  string line, str;
  int nCol = nCols(input_name);
  int nRow = nRows(input_name);
  double sos[nCol],m[nCol],f;
  for (int x = 0; x < nCol; x++){
    sos[x] = 0;
    m[x] = 0;
  }
  compute_mean(input_name,m);
  ifstream file(input_name);
  while(getline(file,line) && !line.empty()) {
    istringstream issline(line);
    int cid = 0;
    while(getline(issline,str,',')){
      f = stof(str,&sz);
      sos[cid]+= powf(f-m[cid],2);
    //  cout << sos[cid] << " ";
      cid++;
    }
    issline.clear();
  //  cout << endl;
  }
  file.close();
  for (int X = 0; X < nCol; X++){
    sd[X] = sqrt(sos[X]/(nRow-1));
  //  cout << sd[X] << " ";
  }
  //cout << endl;
  //cout << "nRow was: " << nRow <<". nCol was " << nCol << "." <<endl;
  return sd;
}
double* compute_CI(char* input_name,double CI[], double alpha)
{
  int nCol = nCols(input_name);
  int nRow = nRows(input_name);
  double mean[nCol], sd[nCol];
  compute_mean(input_name,mean);
  compute_standard_deviation(input_name,sd);
  normal_distribution<> dist;
  double Z = quantile(dist, 1-alpha / 2);
  //cout << "Z score is:" << T << endl;
  double w = 0;
  for (int X =0; X < nCol; X++){
    w = Z * sd[X] / sqrt(double(nRow));
    CI[0 + 2*X] = mean[X] - w;
    CI[1 + 2*X] = mean[X] + w;
    //cout << "{"<< CI[0 + 2*X] << "," <<  CI[1 + 2*X] << "}" << endl;
  }
  return CI;
}
void print_mean(char* mean_file,double mean[],int nCol)
{
ofstream OSTRM(mean_file);
OSTRM.precision(double_DECIMALS);
OSTRM.setf(ios::fixed);
string sep = "";
for (int X = 0;X < nCol; X++) {
  OSTRM << sep << mean[X];
  sep=",";
}
OSTRM << endl;
OSTRM.close();
}
void print_standard_deviation(char* sd_file,double sd[],int nCol)
{
  ofstream OSTRM(sd_file);
  OSTRM.precision(double_DECIMALS);
  OSTRM.setf(ios::fixed);
  string sep = "";
  for (int X = 0;X < nCol; X++) {
    OSTRM << sep << sd[X];
    sep=",";
  }
  OSTRM << endl;
  OSTRM.close();
}
void print_CI(char* CI_file,double CI[], int nCol)
{
  ofstream OSTRM(CI_file);
  OSTRM.precision(double_DECIMALS);
  OSTRM.setf(ios::fixed);
  for (int Y = 0; Y<2; Y++){
    string sep = "";
    for (int X = 0;X < nCol; X++){
      OSTRM << sep << CI[Y+ 2*X];
      sep=",";
    }
    OSTRM << endl;
  }
  OSTRM.close();
}
int main(int argc, char **argv)
{
if (argc != 6) {
  cout << "Usage: ./stats input-file output-mean-file output-sd-file output-CI-file alpha" << endl; //output-mean-file output-sd-file output-CI95-file
  exit(0);
}
int nCol = nCols(argv[1]);
//int nRow = nRows(argv[1]);
double alpha = (double)atoi(argv[5])/100;
//cout << "alpha is: " << alpha << endl;
double mean[nCol], sd[nCol], CI[2 * nCol];
compute_mean(argv[1],mean);
compute_standard_deviation(argv[1],sd);
compute_CI(argv[1],CI,alpha);
print_mean(argv[2],mean,nCol);
print_standard_deviation(argv[3],sd,nCol);
print_CI(argv[4],CI,nCol);
/*
#ifdef USE_OUTPUT_FILE
  ofstream outfile(argv[2]);
  if (!outfile.is_open()) {
    cerr << "Output file isn't open" << endl;
    exit(1);
  }
#endif
#ifdef TIME_PROFILE
  auto start = chrono::high_resolution_clock::now();
#endif
  // OSTRM.precision(double_DECIMALS);
  // OSTRM.setf(ios::fixed);
  // sep = "";
  // for (int id = 0;id < nc; id++) {
  //   OSTRM << sep << Mean[id];
  //   sep=",";
  // }
  // OSTRM << endl;
cout << endl << "There were " << to_string(nLines) << " rows." << endl;
#ifdef TIME_PROFILE
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << endl << "Finding the mean took "
       << duration.count() << " microseconds." << endl;
#endif
*/
return 0;
}
