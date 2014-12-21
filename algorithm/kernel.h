#ifndef ALGORITHMKERNEL_H
#define ALGORITHMKERNEL_H

#include "initialguess.h"

namespace Algorithm
{

class Kernel
{

public:
	typedef enum {
		SUCCESS,
		INFEASIBLE,
		MAXSTEPS,
		CHOLFAILED,
		GAUSSFAILED,
		SMALLSTEPSIZE,
		FAILURE
	} T_ERR;

	Kernel(int N, double hstep, double t0, double * par, double eps, int maxsteps, int maxinfcounter);
	~Kernel();

    void initModule(InitialGuess * init, unsigned int stid);
    bool innerLoopEndModule(unsigned int stid, Kernel::T_ERR & err);
    bool outerLoopEndModule(unsigned int stid);
    void prepareSigRhsModule(unsigned int stid);
    T_ERR innerBlockUModule(unsigned int stid, int i);
    T_ERR innerBlockXModule(unsigned int stid, int i);
    T_ERR innerBlockMergeModule(unsigned int stid, int i);
    T_ERR mergeBlocksModule(unsigned int stid);
    T_ERR coreModule(unsigned int stid); 

	T_ERR run(InitialGuess * init, double * x);

	void getWarmstartData(double * z, double * y);
	bool getBFGSData(double * s_b, double * y_b);

	int xSize() const
	{
		return xsize;
	}

	int gSize() const
	{
		return gsize;
	}

	int hSize() const
	{
		return hsize;
	}

	int getOuterSteps() const
	{
		return outersteps[0];
	}

	int getSteps() const
	{
		return steps[0];
	}

	double getEps() const
	{
		return currenteps[0];
	}

	double getMu() const
	{
		return mu[0];
	}


	/* all Problems */
	int N;
	double hstep;
    double t0;

    void(*func_rhs)(double * t, double * x, double * u, double * p, double * result);
    void(*func_rhsx)(double * t, double * x, double * u, double * p, double * result);
    void(*func_rhsu)(double * t, double * x, double * u, double * p, double * result);
    void(*func_q_u)(int * i, int * N, double * u, double * p, double * res);
    void(*func_dq_u)(int * i, int * N, double * u, double * p, double * res);
    void(*func_q_x)(int * i, int * N, double * x, double * p, double * res);
    void(*func_dq_x)(int * i, int * N, double * x, double * p, double * res);
    void(*func_r_u)(int * i, int * N, double * u, double * p, double * res);
    void(*func_dr_u)(int * i, int * N, double * u, double * p, double * res);
    void(*func_r_x)(int * i, int * N, double * x, double * p, double * res);
    void(*func_dr_x)(int * i, int * N, double * x, double * p, double * res);

	int dimx;
	int dimu;
	int dimk;
	int stages;
	double * par;

	double * rkA;
	double * rkb;
	double * rkc;
    double * lb_u;
    int * nlb_u;
    int * plb_u;
    int * ptlb_u;
    double * ub_u;
    int * nub_u;
    int * pub_u;
    int * ptub_u;
    int * nlb_x;
    int * plb_x;
    int * ptlb_x;
    int * nub_x;
    int * pub_x;
    int * ptub_x;
    int * nq_u;
    int * nr_u;
    int * nq_x;
    int * nr_x;
    int total_nq_u;
    int total_nq_x;
    int total_nr_u;
    int total_nr_x;

    int * val_idx_lb_u;
    int * val_idx_ub_u;
    int * val_idx_lb_x;
    int * val_idx_ub_x;
    int * idx_dr_u;
    int * idx_dr_x;
    int * idx_dq_u;
    int * idx_dq_x;

    int xsize;
    int gsize;
    int hsize;
    int nrhs;
    int total_sigu;
    int total_sigx;
    int * idx_sigu;
    int * idx_sigx;

    int maxsteps;
    int maxinfcounter;
    double eps;

    int lwork;
    int * P;
    int * blocksize;
    int * blockidx;
    int * idx_ub_u;
    int * idx_lb_u;
    int * idx_r_u;
    int * idx_ub_x;
    int * idx_lb_x;
    int * idx_r_x;
    int * idx_hk;
    int * idx_hx;
    int * idx_q_u;
    int * idx_q_x;
    int * idx_u;
    int * idx_x;
    int * idx_k;

    int * idx_y24;
    int * idx_y6810;

    int * merge_sbs;
    int * merge_ss3;
    int * merge_ss6;
    int * merge_sdiag;
    int * merge_nbs;
    int * merge_ns3;
    int * merge_ns6;

    int ** sigucholscheduler;
    int ** sigxcholscheduler;
    int ** mergecholscheduler;

    /* single Problem */
    double ** lb_x;
    double ** ub_x;

    int size_dfu;
    int size_dfx;
    int size_dqu;
    int size_dqx;
    int size_sigu;
    int size_sigx;
    int size_y13;
    int size_y24;
    int size_y579;
    int size_y6810;

    double * xi;
    double ** g;
    double ** h;
    double * dfu;
    double * dfx;
    double ** dfxA;
    double * dqu;
    double * dqx;
    double ** dru;
    double ** drx;
    double ** vec_x;
    double ** vec_s;
    double ** vec_x_tmp;
    double ** vec_s_tmp;
    double ** vec_y;
    double ** vec_z;
    double ** p_s;
    double ** p_z;
    double ** s_b;
    double ** y_b;
    double ** s_b_tmp;
    double ** y_b_tmp;
    double ** sig;
    double * sigx;
    double * sigu;
    int * Bcomputed;
    double * mu;
    double * minustau;
    double * currenteps;
    int * steps;
    int * outersteps;
    int * infcounter;
    double * cx1;
    double * cx2;
    double * phi0;
    double * phia;
    double * dphi0;
    double * alpha_s;
    double * alpha_z;
    int * returncode;
    int * grididx;

    double ** lgsb;
    double * work;
    double * y13;
    double * y24;
    double * y579;
    double * y6810;
    double ** Q1;
    double ** Q2;
    double ** R;
    double ** DN;
    double ** cN;

protected:

    void evrest(int stid, int dograd);
    
    /* chol-scheduler */
    typedef struct { int data[3]; } T_CONV_3;
	typedef struct { int data[6]; } T_CONV_6;

	bool convertIndex(int & row, int & col);
	void generateCholScheduler();
    void generateSigCholScheduler();

};


}

#endif
