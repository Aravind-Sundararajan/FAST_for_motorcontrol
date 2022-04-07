#ifndef BASE_H_
#define BASE_H_
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <gmp.h>
#include <float.h>
#include <math.h>
#include <limits>
#include <numeric>
#include <chrono>
#include <random>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <eigen3/Eigen/Dense>
#include "ball_walker.hpp"
#include "dikin_walker.hpp"
#include "vaidya_walker.hpp"
#include "john_walker.hpp"
#include "public_interface.hpp"
#include "fusion.h"
extern "C"
{
	#include "lrsdriver.h"
	#include "lrslib.h"
}
#define OOB (std::size_t)-1
typedef std::numeric_limits< double > dbl;
#define BIGM dbl::max()
//#define HARN (P.H.size_y+P.H.size_x)*pow(P.H.size_x-1,3)// HARN^2 is the number of steps we should run HAR for in window.h
#define TOL 0.000000000000001 // eps value for  checking equality with zero
#define FUDGE 0.001 // how much to fudge the equality constrain?
#define MAXITER 100000
#define MAXCOL 1000 //LRS
#define FLOAT_DECIMALS 7
#define INT_DECIMALS   7
#define USE_OUTPUT_FILE //write to file in driver programs //probably dont comment this
//#define TIME_PROFILE // uncomment to time program
//#define DEBINFO//uncomment to show debug info
#endif
