#ifndef HEADERFILE_HAR_
#define HEADERFILE_HAR_
#include "base.h"
#include "util.h"
#include "kvp.h"
#include "vec.h"
#include "mat.h"
#include "polytope.h"
#include "optim.h"
#include "redund.h"
using namespace std;
using namespace mosek::fusion;
using namespace monty;
template <typename T>
class hitandrun
{
public:
	polytope<T> P;
	std::size_t n_samples;
	std::size_t thin;
	vec<T> direction;
	vec<T> current;
	std::mt19937 mt;
	std::uniform_real_distribution<T> udist; //we need one uniform [0,1] dist
	std::normal_distribution<T> ndist; // we need one normal  [-inf,inf] dist
	hitandrun()
	{
		polytope<T> P;
		n_samples=0;
		thin=0;
		current   = vec<T>(1);
		direction   = vec<T>(1);
	};
	hitandrun(const hitandrun<T>& HAR)
	{
		P = HAR.P;
		n_samples = HAR.n_samples;
		thin = HAR.thin;
		current = HAR.current;
		direction = HAR.direction;
	};
	hitandrun<T>& operator=(const hitandrun<T>& HAR)
	{
		P = HAR.P;
		n_samples = HAR.n_samples;
		thin = HAR.thin;
		current = HAR.current;
		direction = HAR.direction;
		return (*this);
	};
	hitandrun(const polytope<T>& P_in, const vec<T>& starting_point_in, std::size_t n_samples_in, std::size_t thin_in)
	{
		P = P_in;
		current   = starting_point_in;
		n_samples = n_samples_in;
		thin = thin_in;
		direction   = vec<T>(current.size);
	};
	~hitandrun(){};
	// init the random device and create instance ndist and dist
	//we may be able to use the normal dist mapped using
	void init()
	{
		std::size_t t = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		mt.seed(t);
		udist = std::uniform_real_distribution<T>(0,1);
		ndist = std::normal_distribution<T>(0,1); // u = 0, sigma = 1!
	};
	//investigations into rcond found
	//the function e^x / (1 + e^x) has a nice property of mapping the reals to the domain (0,1)
	//so we can sample any arbitrary pdf (-inf, inf) ->  (0,1)
	// T dist_mapped(auto distr,std::mt19937 mt){
	// 	T v = distr(mt);
	// 	return exp(v) / (1 + exp(v));
	// }
	//main sampling method, writes to file
	void run(const std::string& out_name)
	{
		init();
		ofstream OSTRM(out_name);
		if (!OSTRM.is_open()) {
			cerr << "Output file isn't open" << endl;
			exit(1);
		}
		for (std::size_t t = 0; t < n_samples*thin; t++){
			if (t % thin == 0){
				print_current(OSTRM);
			}
			set_direction();
			P.find_auxilliar_points();
			P.compute_lambdas(current,direction);
			step(direction);
		}
		OSTRM.close();
	};
	//sampling method just iters, no return
	void run()
	{
		init();
		for (std::size_t t = 0; t < n_samples*thin; t++){
			while (step(direction) == 1);
			// print_direction();
			// print_current();
		}
	};
	mat<T> run(const mat<T> out)
	{
		init();
		for (std::size_t t = 0; t < n_samples*thin; t++){
			while (step(direction) == 1);
			if (t % thin == 0) out = out.concat(current);
		}
		return out;
	};

	int step(const vec<T>& obj)
	{
		set_direction();
		P.find_auxilliar_points();
		P.compute_lambdas(current,direction);
		T lampos =  FLT_MAX;//lambda
		T lamneg =  -FLT_MAX;
		for (std::size_t Y =0; Y < P.size_y; Y++){
			if (P.lambdas(Y) != NAN){
				if (P.lambdas(Y) > 0){
					if (P.lambdas(Y) < lampos){
						lampos = P.lambdas(Y);
					}
				}
				if (P.lambdas(Y) < 0){
					if (P.lambdas(Y) > lamneg){
						lamneg = P.lambdas(Y);
					}
				}
			}
		}
		if ((lampos == FLT_MAX) && (lamneg == -FLT_MAX)){
			//cout << "dir doesn't intersect any halfspaces." << endl;
			return 1;
		}
		if (lampos == FLT_MAX){
			lampos =0;
		}
		if (lamneg == -FLT_MAX){
			lamneg =0;
		}
		std::uniform_real_distribution<T> ldist(lamneg,lampos);
		T lam = ldist(mt);//T lam = lamneg + (lampos - lamneg)*udist(mt);
		current+=obj*lam;
		return 0;
	};

int chebyshevCenter()
{
	//chebyshevs center
	vec<T> obj(P.H.size_x+1); obj.set_all(0); obj(obj.size-1) = -1;
	mat<T> Aim = P.A.fun(pow2);
	// cout << "aim" << endl;
	// Aim.print();
	vec<T> Ai = Aim.sum(0);
	Ai = Ai.fun(powroot);
	// cout << "ai" << endl;
	// Ai.print();
	mat<T> aug = ((P.H.transpose()).concat(-Ai)).transpose();
	vec<T> opt(P.H.size_x);
	opt = aug.lpsolve(obj);
	current = opt.trunc(1,opt.size-2); //remove the aux var and the lp opt
	// cout << "current" << endl;
	// current.print();
	return 0;
}

int analyticalCenter()
{//conic optimization by the logarithmic potential function to find analytic center of polytope P
	vector< vector<T> > Av(P.A.size_y, vector<T> (P.A.size_x, 0));
	vector< T > bv(P.A.size_y,0);
	for (std::size_t Y =0; Y< P.A.size_y; Y++){
		bv[Y] = P.b(Y);
		for (std::size_t X =0; X< P.A.size_x; X++){
			Av[Y][X] = P.A(Y,X);
		}
	}
	Model::t M = new Model("test");
	auto _M = finally([&]() { M->dispose();});
	//M->setLogHandler([ = ](const std::string & msg) { std::cout << msg << std::flush; } );
	Variable::t x = M->variable("x", P.A.size_x);
	auto A = monty::new_array_ptr<T>(Av);
	auto b = monty::new_array_ptr<T>(bv);
	int k = A->size(0);
	auto u = M->variable(k);
	M->constraint(Expr::hstack(Expr::sub(b,Expr::mul(A, x)),Expr::constTerm(k, 1.0),u), Domain::inPExpCone());
	M->objective("obj", ObjectiveSense::Maximize, Expr::sum(u));
	try {
	M->solve();
	M->acceptedSolutionStatus(AccSolutionStatus::Optimal);
	auto sol = *(x->level());
	for (int X =0; X< P.A.size_x; X++){
		current(X) = sol[X];
	}
	return 0;
	}
	catch (const SolutionError& e){
    // The solution with at least the expected status was not available.
    // We try to diagnoze why.
    std::cout << "Requested solution was not available.\n";
    auto prosta = M->getProblemStatus();
    switch(prosta)
    {
      case ProblemStatus::DualInfeasible:
        std::cout << "Dual infeasibility certificate found.\n";
        break;

      case ProblemStatus::PrimalInfeasible:
        std::cout << "Primal infeasibility certificate found.\n";
        break;

      case ProblemStatus::Unknown:
        // The solutions status is unknown. The termination code
        // indicates why the optimizer terminated prematurely.
        std::cout << "The solution status is unknown.\n";
        char symname[MSK_MAX_STR_LEN];
        char desc[MSK_MAX_STR_LEN];
        MSK_getcodedesc((MSKrescodee)(M->getSolverIntInfo("optimizeResponse")), symname, desc);
        std::cout << "  Termination code: " << symname << " " << desc << "\n";
        break;

      default:
        std::cout << "Another unexpected problem status: " << prosta << "\n";
    }
  }
  catch (const std::exception& e){
    std::cerr << "Unexpected error: " << e.what() << "\n";
  }


	// for (std::size_t Y =0; Y< P.A.size_y; Y++){
	// 	cout << bv[Y] << " ";
	// 	for (std::size_t X =0; X< P.A.size_x; X++){
	// 		cout << -Av[Y][X] << " ";
	// 	}
	// 	cout << endl;
	// }
	M->dispose();
	  return 1;
}

int analyticalCenterNewtons()
{
	mat<T> y2(P.b.size,1);
	// cout << "y2" << endl;
	// y2.print();
	mat<T> y(P.b.size,1);
	// cout << "y" << endl;
	// y.print();
	for (std::size_t t = 0; t < 100;t++){
		for (std::size_t Y = 0; Y < P.b.size; Y++){
			y(Y,0) = P.b(Y) - P.A.row(Y) % current;
			y2(Y,0) = pow(y(Y,0),2);
		}
		mat<T> y2d = make_diag<T>(y2.col(0));
		// cout << "y2d" << endl;
		// y2d.print();
		mat<T> hess = P.A.transpose() % y2d % P.A;
		// cout << "hess" << endl;
		// hess.print();
		mat<T> hess_inv = hess.inverse();
		// cout << "hessinv" << endl;
		// hess_inv.print();
		mat<T> dk = -hess_inv % P.A.transpose() % y;
		// cout << "dk" << endl;
		// dk.print();
		for (std::size_t X = 0; X < current.size; X++){
			current(X) = current(X) + .9*dk(X,0);
		}
		if (abs(dk.frobenius()) < TOL){
			cout << "current: " << t  << ":" << dk.frobenius()<< endl;
			current.print();
			return 0;
		}
	}
	return 1;
}

	void set_direction()
	{
		for (std::size_t X = 0; X< direction.size; X++){
			direction(X) = ndist(mt); //normal dist domain [-inf, inf]
		}
		vec<T> sos = direction.fun(pow2);
		T norm2 = sos.sum();
		T norm = sqrt(norm2);
		direction/=norm;
	};
	void print_current(std::ofstream& OSTRM)
	{
		std::size_t X;
		string sep;
		sep = "";
		for (X = 0; X< P.size_x; X++){
			OSTRM << sep << current(X);
			sep=",";
		}
		OSTRM << endl;
	};

		void print_current(){
			cout<< "current:[" << current.size << "]" << endl;
			current.print();
		};
		void print_direction(){
			cout<< "direction:[" << direction.size<< "]" << endl;
			direction.print();
		};
		void debinfo(std::size_t t){
			cout << "sample:" << t << endl;
			print_current();
			print_direction();
			P.print_aux();
			P.print_lambdas();
		};
		void dump(){
			cout << "dumping hit and run:" << endl;
			P.dump();
			print_current();
			print_direction();
		};
#ifdef DEBINFO
	#endif
};
#endif
