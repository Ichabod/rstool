#include <list>
#include <cstring>
#include <cmath>
#include <iostream>
#include <set>

#include "kernel.h"
#include "initialguess.h"
#include "../deosapp/system.h"
#include "../scheduler/storage.h"
#include "linalg.h"
#include "fortranheader.h"
#include "../optchol/fortranheader.h"
#include "../optchol/scheduler.h"

#define CONDFREE(p) { if (p!=0) delete[] p; }

using namespace std;

namespace Algorithm
{

template<typename T>
T ** storage_alloc(int size)
{
	T ** res = new T*[STORAGE_SIZE];
	for (int i=0; i<STORAGE_SIZE; i++)
	{
		res[i] = new T[size];
	}
	return res;
}

template<typename T>
T * linear_storage_alloc(int size)
{
	return new T[size*STORAGE_SIZE];
}

template<typename T>
void storage_free(T ** ptr)
{
	for (int i=0; i<STORAGE_SIZE; i++)
	{
		delete[] ptr[i];
	}
	delete[] ptr;
}

template<typename T>
void linear_storage_free(T * ptr)
{
	delete[] ptr;
}

Kernel::Kernel(int N, double hstep, double t0, double * par, double eps, int maxsteps, int maxinfcounter)
{
	//Konstanten
	this->par = par;
	this->N=N;
	this->hstep=hstep;
	this->t0=t0;
	this->eps=eps;
	this->maxsteps = maxsteps;
	this->maxinfcounter = maxinfcounter;

	// System;
	func_rhs = rhs_;
	func_rhsx = rhsx_;
	func_rhsu = rhsu_;
	func_q_u = q_u_;
	func_dq_u = dq_u_;
	func_q_x = q_x_;
	func_dq_x = dq_x_;
	func_r_u = r_u_;
	func_dr_u = dr_u_;
	func_r_x = r_x_;
	func_dr_x = dr_x_;

	readnx_(&(dimx));
	readnu_(&(dimu));
	rdrks_(&(stages));
	dimk=dimx*stages;
	rkA=new double [stages*stages];
	rkb=new double [stages];
	rkc=new double [stages];
	readrk_(rkA,rkb,rkc);

	nlb_u=new int[N];
	nub_u=new int[N];
	nq_u=new int[N];
	nr_u=new int[N];
	val_idx_lb_u=new int[N];
	val_idx_ub_u=new int[N];
	idx_dq_u=new int[N];
	idx_dr_u=new int[N];
	nbox_u_(&N, par, nlb_u, nub_u);
	nq_u_(&N, par, nq_u);
	nr_u_(&N, par, nr_u);
	int total_nlb_u=0;
	int total_nub_u=0;
	total_nq_u=0;
	total_nr_u=0;
	val_idx_lb_u[0] = 1;
	val_idx_ub_u[0] = 1;
	idx_dq_u[0]=0;
	idx_dr_u[0]=0;
	for (int i=1; i<=N; i++)
	{
		total_nlb_u += nlb_u[i-1];
		total_nub_u += nub_u[i-1];
		total_nq_u += nq_u[i-1];
		total_nr_u += nr_u[i-1];
		if (i<N)
		{
			val_idx_lb_u[i]=val_idx_lb_u[i-1]+nlb_u[i-1];
			val_idx_ub_u[i]=val_idx_ub_u[i-1]+nub_u[i-1];
			idx_dq_u[i]=idx_dq_u[i-1] + nq_u[i-1]*dimu;
			idx_dr_u[i]=idx_dr_u[i-1] + nr_u[i-1]*dimu;
		}
	}
	plb_u=new int[total_nlb_u];
	pub_u=new int[total_nub_u];
	ptlb_u=new int[dimu*N];
	ptub_u=new int[dimu*N];
	lb_u=new double [total_nlb_u];
	ub_u=new double [total_nub_u];
	for (int i=1; i<=N; i++)
	{
		pbox_u_(&i, &N, par, &(plb_u[val_idx_lb_u[i-1]-1]), &(pub_u[val_idx_ub_u[i-1]-1]));

		for (int j=0; j<dimu; j++)
		{
			ptlb_u[(i-1)*dimu + j]=-1;
			ptub_u[(i-1)*dimu + j]=-1;

			for (int k=0; k<nlb_u[i-1]; k++)
			{
				if (plb_u[val_idx_lb_u[i-1]-1+k] == j+1)
				{
					ptlb_u[(i-1)*dimu + j]=k+1;
					break;
				}
			}
			for (int k=0; k<nub_u[i-1]; k++)
			{
				if (pub_u[val_idx_ub_u[i-1]-1+k] == j+1)
				{
					ptub_u[(i-1)*dimu + j]=k+1;
					break;
				}
			}
		}

		box_u_(&i, &N, par, &(lb_u[val_idx_lb_u[i-1]-1]), &(ub_u[val_idx_ub_u[i-1]-1]));
	}

	nlb_x=new int[(N+1)];
	nub_x=new int[(N+1)];
	nq_x=new int[(N+1)];
	nr_x=new int[(N+1)];
	val_idx_lb_x=new int[(N+1)];
	val_idx_ub_x=new int[(N+1)];
	idx_dq_x=new int[(N+1)];
	idx_dr_x=new int[(N+1)];
	nbox_x_(&N, par, nlb_x, nub_x);
	nq_x_(&N, par, nq_x);
	nr_x_(&N, par, nr_x);
	int total_nlb_x=0;
	int total_nub_x=0;
	total_nq_x=0;
	total_nr_x=0;
	val_idx_lb_x[0] = 1;
	val_idx_ub_x[0] = 1;
	idx_dq_x[0]=0;
	idx_dr_x[0]=0;
	for (int i=1; i<=N+1; i++)
	{
		total_nlb_x += nlb_x[i-1];
		total_nub_x += nub_x[i-1];
		total_nq_x += nq_x[i-1];
		total_nr_x += nr_x[i-1];
		if (i<N+1)
		{
			val_idx_lb_x[i]=val_idx_lb_x[i-1]+nlb_x[i-1];
			val_idx_ub_x[i]=val_idx_ub_x[i-1]+nub_x[i-1];
			idx_dq_x[i]=idx_dq_x[i-1] + nq_x[i-1]*dimx;
			idx_dr_x[i]=idx_dr_x[i-1] + nr_x[i-1]*dimx;
		}
	}
	plb_x=new int[total_nlb_x];
	pub_x=new int[total_nub_x];
	ptlb_x=new int[dimx*(N+1)];
	ptub_x=new int[dimx*(N+1)];
	lb_x=storage_alloc<double>(total_nlb_x);
	ub_x=storage_alloc<double>(total_nub_x);
	for (int i=1; i<=N+1; i++)
	{
		pbox_x_(&i, &N, par, &(plb_x[val_idx_lb_x[i-1]-1]), &(pub_x[val_idx_ub_x[i-1]-1]));

		for (int j=0; j<dimx; j++)
		{
			ptlb_x[(i-1)*dimx + j]=-1;
			ptub_x[(i-1)*dimx + j]=-1;

			for (int k=0; k<nlb_x[i-1]; k++)
			{
				if (plb_x[val_idx_lb_x[i-1]-1+k] == j+1)
				{
					ptlb_x[(i-1)*dimx + j]=k+1;
					break;
				}
			}
			for (int k=0; k<nub_x[i-1]; k++)
			{
				if (pub_x[val_idx_ub_x[i-1]-1+k] == j+1)
				{
					ptub_x[(i-1)*dimx + j]=k+1;
					break;
				}
			}
		}

		box_x_(&i, &N, par, &(lb_x[0][val_idx_lb_x[i-1]-1]), &(ub_x[0][val_idx_ub_x[i-1]-1]));
	}

	xsize = N*dimu+(N+1)*dimx + N*stages*dimx;
	gsize = total_nub_u + total_nlb_u + total_nr_u + total_nub_x + total_nlb_x + total_nr_x;
	hsize = (stages+1)*dimx*N + total_nq_u + total_nq_x;
	nrhs = xsize+hsize;

	g=storage_alloc<double>(gsize);
	h=storage_alloc<double>(hsize);
	size_dfu = N*stages*dimu*dimx;
	dfu=linear_storage_alloc<double>(size_dfu);
	size_dfx = N*stages*dimx*dimx;
	dfx=linear_storage_alloc<double>(size_dfx);
	dfxA=storage_alloc<double>(N*stages*dimx*stages*dimx);
	size_dqu=total_nq_u*dimu;
	dqu=linear_storage_alloc<double>(size_dqu);
	size_dqx=total_nq_x*dimx;
	dqx=linear_storage_alloc<double>(size_dqx);
	dru=storage_alloc<double>(total_nr_u*dimu);
	drx=storage_alloc<double>(total_nr_x*dimx);

	vec_x=storage_alloc<double>(xsize);
	vec_s=storage_alloc<double>(gsize);
	vec_x_tmp=storage_alloc<double>(xsize);
	vec_s_tmp=storage_alloc<double>(gsize);
	vec_y=storage_alloc<double>(hsize);
	vec_z=storage_alloc<double>(gsize);
	p_s=storage_alloc<double>(gsize);
	p_z=storage_alloc<double>(gsize);
	s_b=storage_alloc<double>(xsize);
	y_b=storage_alloc<double>(xsize);
	s_b_tmp=storage_alloc<double>(xsize);
	y_b_tmp=storage_alloc<double>(xsize);
	idx_sigu=new int [N];
	idx_sigx=new int [(N+1)];
	idx_sigu[0]=1;
	idx_sigx[0]=1;
	total_sigu=0;
	total_sigx=0;
	for (int i=1; i<N; i++)
	{
		if (nr_u[i-1]==0)
		{
			// Nur die Diagonale
			idx_sigu[i] = idx_sigu[i-1] + dimu;
			total_sigu+= dimu;
		}
		else
		{
			// ganze Matrix
			idx_sigu[i] = idx_sigu[i-1] + dimu*dimu;
			total_sigu+= dimu*dimu;
		}
	}
	total_sigu += (nr_u[N-1]==0?dimu:dimu*dimu);
	for (int i=1; i<N+1; i++)
	{
		if (nr_x[i-1]==0)
		{
			idx_sigx[i] = idx_sigx[i-1] + dimx;
			total_sigx+= dimx;
		}
		else
		{
			idx_sigx[i] = idx_sigx[i-1] + dimx*dimx;
			total_sigx+= dimx*dimx;
		}
	}
	total_sigx += (nr_x[N]==0?dimx:dimx*dimx);

	size_sigu=total_sigu;
	sigu=linear_storage_alloc<double>(total_sigu);
	size_sigx=total_sigx;
	sigx=linear_storage_alloc<double>(total_sigx);
	sig=storage_alloc<double>(gsize);

	lgsb=storage_alloc<double>((xsize+hsize)*3);
	P=new int [(xsize+hsize)];
	blocksize=new int [(N+2)];
	blockidx=new int [(N+2)];
	idx_ub_u=new int [N];
    idx_lb_u=new int [N];
    idx_r_u=new int [N];
    idx_ub_x=new int [(N+1)];
    idx_lb_x=new int [(N+1)];
    idx_r_x=new int [(N+1)];
    idx_hk=new int [stages*N];
    idx_hx=new int [N];
    idx_q_u=new int [N];
    idx_q_x=new int [(N+1)];
    idx_u=new int [N];
    idx_x=new int [(N+1)];
    idx_k=new int [stages*N];

    lwork = (xsize+hsize)*3;
    work = linear_storage_alloc<double>(lwork);

    size_y13 = dimu*N*(dimx*stages+3);
	y13 = linear_storage_alloc<double>(size_y13);
	size_y24 = total_nq_u*(dimx*stages+3);
	y24 = linear_storage_alloc<double>(size_y24);
	size_y579 = dimx*N*(dimx*stages+dimx+3);
	y579 = linear_storage_alloc<double>(size_y579);
	size_y6810 = (total_nq_x-nq_x[N])*(dimx*stages+dimx+3);
	y6810 = linear_storage_alloc<double>(size_y6810);

	idx_y24=new int [N*sizeof(int)];
	idx_y6810=new int [N*sizeof(int)];
	idx_y24[0] = 1;
	idx_y6810[0] = 1;
	for (int i=1; i<N; i++)
	{
		idx_y24[i] = idx_y24[i-1] + nq_u[i-1]*(dimx*stages+3);
		idx_y6810[i] = idx_y6810[i-1] + nq_x[i-1]*(dimx*stages+dimx+3);
	}



	int index;

	index=0;
    for (int i=0; i<N; i++)
    {
        idx_ub_u[i] = index;
        index += nub_u[i];
        idx_lb_u[i] = index;
        index += nlb_u[i];
        idx_r_u[i] = index;
        index += nr_u[i];
        
        idx_u[i] = i*dimu;
    }
    for (int i=0; i<N+1; i++)
    {
        idx_ub_x[i] = index;
        index = index + nub_x[i];
        idx_lb_x[i] = index;
        index = index + nlb_x[i];
        idx_r_x[i] = index;
        index = index + nr_x[i];

        idx_x[i] = N*dimu + i*dimx;
    }

 	index=0;
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<stages; j++)
        {
            idx_hk[i*stages + j] = index;
            index = index+dimx;
            idx_k[i*stages + j] = N*dimu + (N+1)*dimx + i*stages*dimx + j*dimx;
        }
        idx_hx[i] = index;
        index=index+dimx;
    }
    for (int i=0; i<N; i++)
    {
        idx_q_u[i] = index;
        index = index + nq_u[i];
    }
    for (int i=0; i<N+1; i++)
    {
        idx_q_x[i] = index;
        index = index + nq_x[i];
    }

    index=0;
    for (int i=0; i<N; i++)
    {
    	for (int j=0; j<dimu; j++) P[index+j] = i*dimu+j+1;
        index=index+dimu;
        for (int j=0; j<nq_u[i]; j++) P[index+j] = xsize + idx_q_u[i]+j+1;
        index=index+nq_u[i];
        for (int j=0; j<dimx; j++) P[index+j] = N*dimu + i*dimx+j+1;
        index=index+dimx;
        for (int j=0; j<nq_x[i]; j++) P[index+j] = xsize + idx_q_x[i]+j+1;
        index=index+nq_x[i];
        for (int j=0; j<dimx*stages; j++) P[index+j] = N*dimu + (N+1)*dimx + i*dimx*stages+j+1;
        index=index+dimx*stages;
        for (int j=0; j<dimx*stages; j++) P[index+j] = xsize + idx_hk[i*stages]+j+1;
        index=index+dimx*stages;
    }
    for (int j=0; j<dimx; j++) P[index+j] = N*dimu + N*dimx+j+1;
    index=index+dimx;
    for (int j=0; j<nq_x[N]; j++) P[index+j] = xsize + idx_q_x[N]+j+1;
    index=index+nq_x[N];
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<dimx; j++) P[index+j] = xsize + idx_hx[i]+j+1;
        index=index+dimx;
    }

    blockidx[0] = 0;
    for (int i=0; i<N; i++)
    {
        blocksize[i] = dimu+nq_u[i]+dimx+nq_x[i]+dimx*stages+dimx*stages;
        blockidx[i+1] = blockidx[i]+blocksize[i];
    }
    blocksize[N] = dimx+nq_x[N];
    blocksize[N+1] = N*dimx;
    blockidx[N+1] = blockidx[N]+blocksize[N];

	Q1 = storage_alloc<double>((blockidx[N]+blocksize[N])*dimx);
	R = storage_alloc<double>((blockidx[N]+blocksize[N])*3);
	Q2 = storage_alloc<double>((blockidx[N-1]+blocksize[N-1]-blocksize[0])*dimx);

	DN = storage_alloc<double>((N*dimx*dimx+ (N-1)*dimx*dimx));
	cN = storage_alloc<double>(N*dimx*3);

	xi = new double[STORAGE_SIZE];
	Bcomputed = new int[STORAGE_SIZE];
	mu = new double[STORAGE_SIZE];
	minustau = new double[STORAGE_SIZE];
	currenteps = new double[STORAGE_SIZE];
	steps = new int[STORAGE_SIZE];
	outersteps = new int[STORAGE_SIZE];
	infcounter = new int[STORAGE_SIZE];
	cx1 = new double[STORAGE_SIZE];
	cx2 = new double[STORAGE_SIZE];
	phi0 = new double[STORAGE_SIZE];
	phia = new double[STORAGE_SIZE];
	dphi0 = new double[STORAGE_SIZE];
	alpha_s = new double[STORAGE_SIZE];
	alpha_z = new double[STORAGE_SIZE];
	returncode = new int[STORAGE_SIZE];
	grididx = new int[STORAGE_SIZE];

	generateCholScheduler();
	generateSigCholScheduler();

	InitialGuess::xsize = xsize;
	InitialGuess::gsize = gsize;
	InitialGuess::hsize = hsize;

}


Kernel::~Kernel()
{
	CONDFREE(rkA);
	CONDFREE(rkb);
	CONDFREE(rkc);

	CONDFREE(lb_u);
	CONDFREE(nlb_u);
	CONDFREE(plb_u);
	CONDFREE(ptlb_u);
	CONDFREE(ub_u);
	CONDFREE(nub_u);
	CONDFREE(pub_u);
	CONDFREE(ptub_u);
	storage_free<double>(lb_x);
	CONDFREE(nlb_x);
	CONDFREE(plb_x);
	CONDFREE(ptlb_x);
	storage_free<double>(ub_x);
	CONDFREE(nub_x);
	CONDFREE(pub_x);
	CONDFREE(ptub_x);
	CONDFREE(nq_u);
	CONDFREE(nr_u);
	CONDFREE(nq_x);
	CONDFREE(nr_x);
    CONDFREE(idx_dr_u);
    CONDFREE(idx_dr_x);
    CONDFREE(idx_dq_u);
    CONDFREE(idx_dq_x);

	CONDFREE(val_idx_lb_u);
	CONDFREE(val_idx_ub_u);
	CONDFREE(val_idx_lb_x);
	CONDFREE(val_idx_ub_x);

    storage_free<double>(g);
    storage_free<double>(h);
    linear_storage_free<double>(dfu);
    linear_storage_free<double>(dfx);
    storage_free<double>(dfxA);
    linear_storage_free<double>(dqu);
    linear_storage_free<double>(dqx);
    storage_free<double>(dru);
    storage_free<double>(drx);
    storage_free<double>(vec_x);
    storage_free<double>(vec_s);
    storage_free<double>(vec_x_tmp);
    storage_free<double>(vec_s_tmp);
    storage_free<double>(vec_y);
    storage_free<double>(vec_z);
    storage_free<double>(p_s);
    storage_free<double>(p_z);
    storage_free<double>(s_b);
    storage_free<double>(y_b);
    storage_free<double>(s_b_tmp);
    storage_free<double>(y_b_tmp);
    storage_free<double>(sig);
    linear_storage_free<double>(sigx);
    linear_storage_free<double>(sigu);
    CONDFREE(idx_sigx);
    CONDFREE(idx_sigu);

	storage_free<double>(lgsb);
	CONDFREE(P);
	CONDFREE(blocksize);
	CONDFREE(blockidx);
	CONDFREE(idx_ub_u);
    CONDFREE(idx_lb_u);
    CONDFREE(idx_r_u);
    CONDFREE(idx_ub_x);
    CONDFREE(idx_lb_x);
    CONDFREE(idx_r_x);
    CONDFREE(idx_hk);
    CONDFREE(idx_hx);
    CONDFREE(idx_q_u);
    CONDFREE(idx_q_x);
    CONDFREE(idx_u);
    CONDFREE(idx_x);
    CONDFREE(idx_k);
    linear_storage_free<double>(y13);
	linear_storage_free<double>(y24);
	linear_storage_free<double>(y579);
	linear_storage_free<double>(y6810);
	CONDFREE(idx_y24);
	CONDFREE(idx_y6810);
	storage_free<double>(Q1);
	storage_free<double>(Q2);
	storage_free<double>(R);
	storage_free<double>(DN);
	storage_free<double>(cN);

    CONDFREE(merge_sbs);
    CONDFREE(merge_ss3);
    CONDFREE(merge_ss6);
    CONDFREE(merge_sdiag);
    CONDFREE(merge_nbs);
    CONDFREE(merge_ns3);
    CONDFREE(merge_ns6);

    linear_storage_free<double>(work);

	CONDFREE(xi);
	CONDFREE(Bcomputed);
	CONDFREE(mu);
	CONDFREE(minustau);
	CONDFREE(currenteps);
	CONDFREE(steps);
	CONDFREE(outersteps);
	CONDFREE(infcounter);
	CONDFREE(cx1);
	CONDFREE(cx2);
	CONDFREE(phi0);
	CONDFREE(phia);
	CONDFREE(dphi0);
	CONDFREE(alpha_s);
	CONDFREE(alpha_z);
	CONDFREE(returncode);
	CONDFREE(grididx);

	for (int i=0; i<N; i++)
	{
		CONDFREE(sigucholscheduler[i]);
		CONDFREE(mergecholscheduler[i]);
	}
	CONDFREE(sigucholscheduler);
	CONDFREE(mergecholscheduler);

	for (int i=0; i<N+1; i++)
	{
		CONDFREE(sigxcholscheduler[i]);
	}
	CONDFREE(sigxcholscheduler);

}

bool Kernel::convertIndex(int & row, int & col)
{
   int row_block = row/dimx;
   int col_block = col/dimx;

   if (row_block == col_block)
   {
      row %= dimx;
      return true;
   }
   else if (row_block == col_block+1)
   {
      col += N*dimx;
      row %= dimx;
      return true;
   }
   else
   {
      return false;
   }
}

void Kernel::generateCholScheduler()
{
   int cursor;
   int counter;
   int i1,i2,i3,j1,j2,j3;

   merge_ns3 = new int[N*dimx-1];
   merge_ns6 = new int[N*dimx-1];
   merge_nbs = new int[N*dimx-1];
   merge_sdiag = new int[N*dimx*2];

   list<T_CONV_3> conv_ss3;
   list<T_CONV_3> conv_sbs;
   list<T_CONV_6> conv_ss6;

   for (int i=0; i<N*dimx;i++)
   {
//      cout << "Schritt " << i << endl;

      i1 = i;
      j1 = i;
      convertIndex(i1,j1);
//      cout << "  g(" << i << "," << i << ") ===> " ;
//      cout << "  g(" << i1 << "," << j1 << ")" << endl;
      merge_sdiag[i*2] = i1+1;
      merge_sdiag[i*2+1] = j1+1;

      if (i < N*dimx-1)
      {
         counter = 0;
         int bidx = N*dimx-i-1;
         for (int j=0; j<bidx; j++)
         {
            i1 = bidx;
            j1 = j;

            if (convertIndex(i1,j1))
            {
//               cout << "  b(" << j << ") -= b(" << bidx << ")*g(" << bidx << "," << j << ") ===> " ;
//               cout << "  b(" << j << ") -= b(" << bidx << ")*g(" << i1 << "," << j1 << ")" << endl;
               T_CONV_3 data = {{j,i1,j1}};
               conv_sbs.push_back( data );
               counter++;
            }

         }
         merge_nbs[i] = counter;


         counter = 0;
         for (int j=i+1; j<N*dimx; j++)
         {
            i1 = j;
            j1 = i;

            if (convertIndex(i1,j1))
            {
//               cout << "  b(" << j << ") -= b(" << i << ")*g(" << j << "," << i << ") ===> " ;
//               cout << "  b(" << j << ") -= b(" << i << ")*g(" << i1 << "," << j1 << ")" << endl;
               T_CONV_3 data = {{j,i1,j1}};
               conv_ss3.push_back( data );
               counter++;
            }

         }
         merge_ns3[i] = counter;

         counter = 0;
         for (int j=i+1; j < N*dimx; j++)
         {
            for (int k=i+1; k<=j; k++)
            {
               i1 = j;
               j1 = k;
               i2 = j;
               j2 = i;
               i3 = k;
               j3 = i;

               if (convertIndex(i1,j1) && convertIndex(i2,j2) && convertIndex(i3,j3))
               {
//                  cout << "  g(" << i1 << "," << j1 << ") -= g(" << i2 << "," << j2 << ")*g(" << i3 << "," << j3 << ")  ===> ";
//                  cout << "  g(" << i1 << "," << j1 << ") -= g(" << i2 << "," << j2 << ")*g(" << i3 << "," << j3 << ")" << endl;
                  T_CONV_6 data = {{i1,j1,i2,j2,i3,j3}};
                  conv_ss6.push_back( data );
                  counter++;
               }

            }
         }
         merge_ns6[i] = counter;
      }
   }

   merge_sbs = new int[conv_sbs.size()*3];
   cursor = 0;
   for (list<T_CONV_3>::iterator it = conv_sbs.begin(); it != conv_sbs.end(); ++it)
   {
      for (int j=0; j<3; j++)
      {
         merge_sbs[cursor++]=(*it).data[j]+1;
      }
   }

   merge_ss3 = new int[conv_ss3.size()*3];
   cursor = 0;
   for (list<T_CONV_3>::iterator it = conv_ss3.begin(); it != conv_ss3.end(); ++it)
   {
      for (int j=0; j<3; j++)
      {
         merge_ss3[cursor++]=(*it).data[j]+1;
      }
   }

   merge_ss6 = new int[conv_ss6.size()*6];
   cursor = 0;
   for (list<T_CONV_6>::iterator it = conv_ss6.begin(); it != conv_ss6.end(); ++it)
   {
      for (int j=0; j<6; j++)
      {
         merge_ss6[cursor++]=(*it).data[j]+1;
      }
   }
}

void Kernel::generateSigCholScheduler()
{

	int stid=0;

	double tmpxi = xi[stid];

	xi[stid] = 2.345;
	for (int i=0; i<xsize; i++)
	{
		vec_x[stid][i] = 23.456 + (double)i*0.012;
	}
	for (int i=0; i<gsize; i++)
	{
		vec_s[stid][i] = 13.456 + (double)i*0.12;
		vec_z[stid][i] = 1.56 + (double)i*0.124;
	}
	for (int i=0; i<hsize; i++)
	{
		vec_y[stid][i] = 1.917 + (double)i*0.67;
	}

	evrest(stid,1);

	evsig_(xi, vec_s[stid], vec_z[stid], &N, &dimx, &dimu, &gsize,
	nub_u, nlb_u, nub_x, nlb_x,
	pub_u, plb_u, pub_x, plb_x,
	nr_u, nr_x,
	idx_ub_u, idx_lb_u, idx_r_u, idx_ub_x, idx_lb_x, idx_r_x,
	idx_sigu, idx_sigx, sigu+stid*size_sigu, sigx+stid*size_sigx,
	dru[stid], drx[stid]);

	/**********************
	 * evrhs			  *
	 **********************/
	evrhs_(lgsb[stid], &dimx, &dimu, &dimk, &xsize, &hsize, &N,
	nub_u, nlb_u, nub_x, nlb_x,
	pub_u, plb_u, pub_x, plb_x,
	nr_u, nr_x, dru[stid], drx[stid],
	nq_u, nq_x, dqu+stid*size_dqu, dqx+stid*size_dqx,
	g[stid], h[stid], dfu+stid*size_dfu, dfx+stid*size_dfx, dfxA[stid],
	&stages, rkb, &hstep,
	vec_y[stid], vec_z[stid], vec_s[stid], mu+stid,
	y_b[stid], s_b[stid], Bcomputed+stid, xi+stid);

	int mrhs = Bcomputed[stid]==1?3:1;
	int nrhs = xsize+hsize;
	/**********************
	 * permrhs			  *
	 **********************/
	for (int j=0; j<mrhs; j++)
	{
		for (int i=0; i<nrhs; i++)
		{
			work[stid*lwork+j*nrhs+i] = lgsb[stid][j*nrhs+P[i]-1];
		}	
	}


	sigucholscheduler = new int*[N];
	mergecholscheduler = new int*[N];
	sigxcholscheduler = new int*[N+1];

	for (int i=0; i<N; i++)
	{
		if (nr_u[i] > 0)
		{
			sigucholscheduler[i] = OptChol::Scheduler::generate(dimu, sigu+stid*size_sigu+idx_sigu[i]-1);
		} else
		{
			sigucholscheduler[i] = 0;
		}

		if (nr_x[i] > 0)
		{
			sigxcholscheduler[i] = OptChol::Scheduler::generate(dimx, sigx+stid*size_sigx+stid*size_sigx+idx_sigx[i]-1);
		} else
		{
			sigxcholscheduler[i] = 0;
		}

		int err=0;

		if (nq_u[i] > 0)
		{
			int isq = nr_u[i]>0?1:0;
			//TODO: gesamten work ermitteln und reservieren
			double * tempD = new double[nq_u[i]*nq_u[i]];
			double * temp = new double [nq_u[i]*dimu];

			evbur_(&dimu, &dimx, &stages, sigu+stid*size_sigu+idx_sigu[i]-1,
				&isq, dfu+stid*size_dfu+i*dimu*dimx*stages, &nrhs,
					&mrhs, work+stid*lwork, dqu+stid*size_dqu+idx_dq_u[i], &nq_u[i], 
					&blockidx[i], y13+stid*size_y13+i*dimu*(dimx*stages+3),
					y24+stid*size_y24+idx_y24[i]-1, tempD, temp, &err, sigucholscheduler[i]);
			delete[] tempD;
			delete[] temp;
		}
		else
		{
			int isq = nr_u[i]>0?1:0;

			evbu_(&dimu, &dimx, &stages, sigu+stid*size_sigu+idx_sigu[i]-1,
				&isq, dfu+stid*size_dfu+i*dimu*dimx*stages, &nrhs, 
			 	&mrhs, work+stid*lwork, &blockidx[i], y13+stid*size_y13+i*dimu*(dimx*stages+3),
			 	&err, sigucholscheduler[i]);
		}

		if (nq_x[i] > 0 && !err)
		{
			int isq = nr_x[i]>0?1:0;
			int bidx = blockidx[i] + dimu + nq_u[i];
			//TODO: gesamten work ermitteln und reservieren
			double * tempD = new double[nq_x[i]*nq_x[i]];
			double * temp = new double[nq_x[i]*dimx];

			evbxr_(&dimx, &stages, sigx+stid*size_sigx+idx_sigx[i]-1,
				&isq, dfx+stid*size_dfx+i*dimx*dimx*stages, &nrhs,
	     		&mrhs, work+stid*lwork, dqx+stid*size_dqx+idx_dq_x[i], &nq_x[i], 
	     		&bidx, y579+stid*size_y579+i*dimx*(dimx*stages+dimx+3),
	     		y6810+stid*size_y6810+idx_y6810[i]-1, tempD, temp, &err, sigxcholscheduler[i]);
			delete[] tempD;
			delete[] temp;
		}
		else if (!err)
		{
			int isq = nr_x[i]>0?1:0;
			int bidx = blockidx[i] + dimu + nq_u[i];

			evbx_(&dimx, &stages, sigx+stid*size_sigx+idx_sigx[i]-1,
				&isq, dfx+stid*size_dfx+i*dimx*dimx*stages, &nrhs, 
		     	&mrhs, work+stid*lwork, &bidx, y579+stid*size_y579+i*dimx*(dimx*stages + dimx+3),
		     	&err, sigxcholscheduler[i]);
		}

		double * tempD = new double [dimx*stages*dimx*stages];

		if (!err)
		{
			int iQ2=(i>0?1:0);
			//TODO: gesamten work ermitteln und reservieren

			evmrpr_(&dimu, &dimx, &stages, dfu+stid*size_dfu+i*dimu*dimx*stages,
				dfx+stid*size_dfx+i*dimx*dimx*stages,&dfxA[stid][i*dimx*stages*dimx*stages], xi+stid,
				y13+stid*size_y13+i*dimu*(dimx*stages+3),
				y13+stid*size_y13+i*dimu*(dimx*stages+3)+dimu*dimx*stages,
				y579+stid*size_y579+i*dimx*(dimx*stages+ dimx+3),
				y579+stid*size_y579+i*dimx*(dimx*stages+ dimx+3)+dimx*dimx*stages,
				y579+stid*size_y579+i*dimx*(dimx*stages+ dimx+3)+dimx*dimx*stages + dimx*dimx,
				&blockidx[i], &blocksize[i], work+stid*lwork, &nrhs, &mrhs,
				&nq_u[i], &nq_x[i], rkb, &hstep,
				&Q1[stid][blockidx[i]*dimx], (iQ2==0?0:&Q2[stid][(blockidx[i]-blocksize[0])*dimx]), &R[stid][blockidx[i]*3], &err, tempD, &iQ2);
		} else
		{
			memset(tempD, 1, sizeof(double)*dimx*stages*dimx*stages);
		}

		mergecholscheduler[i] = OptChol::Scheduler::generate(dimx*stages, tempD);

		delete[] tempD;
	}

	if (nr_x[N] > 0)
	{
		sigxcholscheduler[N] = OptChol::Scheduler::generate(dimx, sigx+stid*size_sigx+idx_sigx[N]-1);
	} else
	{
		sigxcholscheduler[N] = 0;
	}

	xi[stid] = tmpxi;

}

void Kernel::evrest(int stid, int dograd)
{
	int lwork=dimx+dimx*dimx+dimx*dimu;
	//TODO: gesamte worklänge ermitteln und vorher allokieren
	double * work = new double[lwork];
	int err;

	evrest_(vec_x[stid], par, &t0, &N, &hstep,
     &dimx, &dimu, &xsize, &gsize, &hsize,
     &stages, rkA, rkb, rkc,
     func_rhs, func_rhsx, func_rhsu, 
     func_r_u, func_dr_u, func_r_x, func_dr_x,
     func_q_u, func_dq_u, func_q_x, func_dq_x,
     nub_u, nlb_u, nub_x, nlb_x,
     val_idx_ub_u, val_idx_lb_u, val_idx_ub_x, val_idx_lb_x,
     pub_u, plb_u, pub_x, plb_x,
     ub_u, lb_u, ub_x[stid], lb_x[stid],
     nr_u, nr_x, nq_u, nq_x,
     idx_dr_u, idx_dr_x, idx_dq_u, idx_dq_x,
     idx_ub_u, idx_lb_u, idx_r_u, idx_ub_x, idx_lb_x, idx_r_x,
     idx_q_u, idx_q_x,
     idx_hk, idx_hx,
     g[stid], h[stid], dfu+stid*size_dfu, dfx+stid*size_dfx, dfxA[stid], dqu+stid*size_dqu, dqx+stid*size_dqx, dru[stid], drx[stid],
     work, &lwork, &dograd, &err);

	delete[] work;
}

void Kernel::initModule(InitialGuess * init, unsigned int stid)
{
	/**********************
	 * init				  *
	 **********************/
	outersteps[stid] = 0;
	steps[stid]=0;

	memcpy(vec_x[stid], init->x, sizeof(double)*xsize);

	bool hasstartvalues = (init->z != 0) && (init->y != 0);
	if (hasstartvalues)
	{
		memcpy(vec_z[stid], init->z, gsize*sizeof(double));
		memcpy(vec_y[stid], init->y, hsize*sizeof(double));
	}

	Bcomputed[stid] = (init->s_b != 0) && (init->y_b != 0);
	if (Bcomputed[stid])
	{
		memcpy(s_b[stid], init->s_b, xsize*sizeof(double));
		memcpy(y_b[stid], init->y_b, xsize*sizeof(double));
	}
	int istart=(hasstartvalues==true?0:1);
	xi[stid]=init->xi;
	initip_(&xsize, &gsize, &hsize,
		&istart,
		mu+stid, minustau+stid, currenteps+stid, steps+stid, outersteps+stid, infcounter+stid, 
		vec_s[stid], vec_z[stid], vec_y[stid]);
	if (init->mu!=0.0) mu[stid]=init->mu;

	/**********************
	 * evrest    		  *
	 **********************/
	 evrest(stid, 1);
	 for (int i=0; i<gsize; i++)
	 {
	 	vec_s[stid][i] = max(0.0, g[stid][i])+10.0*eps;
	 }
	 returncode[stid] = SUCCESS;
}

bool Kernel::innerLoopEndModule(unsigned int stid, Kernel::T_ERR & err)
{

	/**********************
	 * ifin 	   		  *
	 **********************/
	int loop = 1;
	int error = 0;

	ifin_(&gsize, &hsize, g[stid], h[stid], currenteps+stid, &eps, infcounter+stid, &maxinfcounter,
     steps+stid, &maxsteps, mu+stid, &error, &loop);

	if (steps[stid] >= maxsteps)
	{
		err = MAXSTEPS;
	} else
	{
		err = ( error==0 ? SUCCESS : INFEASIBLE );	
	}
	if (err != SUCCESS) returncode[stid] = err;

	return loop==1;
}

bool Kernel::outerLoopEndModule(unsigned int stid)
{

	/**********************
	 * ofin 	   		  *
	 **********************/
	int loop = 1;

	ofin_(&gsize, vec_s[stid], vec_z[stid], mu+stid, minustau+stid, currenteps+stid, &eps,
     steps+stid, &maxsteps, outersteps+stid, &loop);

	return loop==1;
}

void Kernel::prepareSigRhsModule(unsigned int stid)
{
	/**********************
	 * evsig			  *
	 **********************/
	evsig_(xi+stid, vec_s[stid], vec_z[stid], &N, &dimx, &dimu, &gsize,
	nub_u, nlb_u, nub_x, nlb_x,
	pub_u, plb_u, pub_x, plb_x,
	nr_u, nr_x,
	idx_ub_u, idx_lb_u, idx_r_u, idx_ub_x, idx_lb_x, idx_r_x,
	idx_sigu, idx_sigx, sigu+stid*size_sigu, sigx+stid*size_sigx,
	dru[stid], drx[stid]);

	/**********************
	 * evrhs			  *
	 **********************/
	evrhs_(lgsb[stid], &dimx, &dimu, &dimk, &xsize, &hsize, &N,
	nub_u, nlb_u, nub_x, nlb_x,
	pub_u, plb_u, pub_x, plb_x,
	nr_u, nr_x, dru[stid], drx[stid],
	nq_u, nq_x, dqu+stid*size_dqu, dqx+stid*size_dqx,
	g[stid], h[stid], dfu+stid*size_dfu, dfx+stid*size_dfx, dfxA[stid],
	&stages, rkb, &hstep,
	vec_y[stid], vec_z[stid], vec_s[stid], mu+stid,
	y_b[stid], s_b[stid], Bcomputed+stid, xi+stid);

	int mrhs = Bcomputed[stid]==1?3:1;
	int nrhs = xsize+hsize;
	/**********************
	 * permrhs			  *
	 **********************/
	for (int j=0; j<mrhs; j++)
	{
		for (int i=0; i<nrhs; i++)
		{
			work[stid*lwork+j*nrhs+i] = lgsb[stid][j*nrhs+P[i]-1];
		}	
	}

}

Kernel::T_ERR Kernel::innerBlockUModule(unsigned int stid, int i)
{
	int mrhs = Bcomputed[stid]==1?3:1;
	int err = 0;

	if (nq_u[i] > 0)
	{
		/**********************
		 * evbur			  *
		 **********************/
		int isq = nr_u[i]>0?1:0;
		//TODO: gesamten work ermitteln und reservieren
		double * tempD = new double[nq_u[i]*nq_u[i]];
		double * temp = new double [nq_u[i]*dimu];

		evbur_(&dimu, &dimx, &stages, sigu+stid*size_sigu+idx_sigu[i]-1,
			&isq, dfu+stid*size_dfu+i*dimu*dimx*stages, &nrhs,
				&mrhs, work+stid*lwork, dqu+stid*size_dqu+idx_dq_u[i], &nq_u[i], 
				&blockidx[i], y13+stid*size_y13+i*dimu*(dimx*stages+3),
				y24+stid*size_y24+idx_y24[i]-1, tempD, temp, &err, sigucholscheduler[i]);
		delete[] tempD;
		delete[] temp;
	}
	else
	{
		/**********************
		 * evbu				  *
		 **********************/
		int isq = nr_u[i]>0?1:0;

		evbu_(&dimu, &dimx, &stages, sigu+stid*size_sigu+idx_sigu[i]-1,
			&isq, dfu+stid*size_dfu+i*dimu*dimx*stages, &nrhs, 
		 	&mrhs, work+stid*lwork, &blockidx[i], y13+stid*size_y13+i*dimu*(dimx*stages+3),
		 	&err, sigucholscheduler[i]);
	}
	if(err!=0)
	{
		return CHOLFAILED;
	}
	
	return SUCCESS;

}

Kernel::T_ERR Kernel::innerBlockXModule(unsigned int stid, int i)
{
	int mrhs = Bcomputed[stid]==1?3:1;
	int err = 0;

	if (nq_x[i] > 0)
	{
		/**********************
		 * evbxr			  *
		 **********************/
		int isq = nr_x[i]>0?1:0;
		int bidx = blockidx[i] + dimu + nq_u[i];
		//TODO: gesamten work ermitteln und reservieren
		double * tempD = new double[nq_x[i]*nq_x[i]];
		double * temp = new double[nq_x[i]*dimx];

		evbxr_(&dimx, &stages, sigx+stid*size_sigx+idx_sigx[i]-1,
			&isq, dfx+stid*size_dfx+i*dimx*dimx*stages, &nrhs,
     		&mrhs, work+stid*lwork, dqx+stid*size_dqx+idx_dq_x[i], &nq_x[i], 
     		&bidx, y579+stid*size_y579+i*dimx*(dimx*stages+dimx+3),
     		y6810+stid*size_y6810+idx_y6810[i]-1, tempD, temp, &err, sigxcholscheduler[i]);
		delete[] tempD;
		delete[] temp;
	}
	else
	{
		/**********************
		 * evbx				  *
		 **********************/
		int isq = nr_x[i]>0?1:0;
		int bidx = blockidx[i] + dimu + nq_u[i];

		evbx_(&dimx, &stages, sigx+stid*size_sigx+idx_sigx[i]-1,
			&isq, dfx+stid*size_dfx+i*dimx*dimx*stages, &nrhs, 
	     	&mrhs, work+stid*lwork, &bidx, y579+stid*size_y579+i*dimx*(dimx*stages + dimx+3),
	     	&err, sigxcholscheduler[i]);
	}
	if(err!=0)
	{
		return CHOLFAILED;
	}

	return SUCCESS;
}

Kernel::T_ERR Kernel::innerBlockMergeModule(unsigned int stid, int i)
{
	int mrhs = Bcomputed[stid]==1?3:1;
	int err = 0;

	int iQ2=(i>0?1:0);
	//TODO: gesamten work ermitteln und reservieren
	double * tempD = new double [dimx*stages*dimx*stages];

	evmrpr_(&dimu, &dimx, &stages, dfu+stid*size_dfu+i*dimu*dimx*stages,
		dfx+stid*size_dfx+i*dimx*dimx*stages, dfxA[stid]+i*dimx*stages*dimx*stages, xi+stid,
		y13+stid*size_y13+i*dimu*(dimx*stages+3),
		y13+stid*size_y13+i*dimu*(dimx*stages+3)+dimu*dimx*stages,
		y579+stid*size_y579+i*dimx*(dimx*stages+ dimx+3),
		y579+stid*size_y579+i*dimx*(dimx*stages+ dimx+3)+dimx*dimx*stages,
		y579+stid*size_y579+i*dimx*(dimx*stages+ dimx+3)+dimx*dimx*stages + dimx*dimx,
		&blockidx[i], &blocksize[i], work+stid*lwork, &nrhs, &mrhs,
		&nq_u[i], &nq_x[i], rkb, &hstep,
		Q1[stid]+blockidx[i]*dimx, (iQ2==0?0:Q2[stid]+(blockidx[i]-blocksize[0])*dimx), R[stid]+blockidx[i]*3, &err, tempD, &iQ2);

	evmrge_(&dimu, &dimx, &stages, dfxA[stid]+i*dimx*stages*dimx*stages,
		xi+stid, y13+stid*size_y13+i*dimu*(dimx*stages+3),
		y24+stid*size_y24+idx_y24[i]-1,
		y13+stid*size_y13+i*dimu*(dimx*stages+3)+dimu*dimx*stages,
		y24+stid*size_y24+idx_y24[i]-1+nq_u[i]*dimx*stages,
		y579+stid*size_y579+i*dimx*(dimx*stages+ dimx+3),
		y6810+stid*size_y6810+idx_y6810[i]-1,
		y579+stid*size_y579+i*dimx*(dimx*stages+ dimx+3)+dimx*dimx*stages,
		y6810+stid*size_y6810+idx_y6810[i]-1 + nq_x[i]*dimx*stages,
		y579+stid*size_y579+i*dimx*(dimx*stages+ dimx+3)+dimx*dimx*stages + dimx*dimx,
		y6810+stid*size_y6810+idx_y6810[i]-1 + nq_x[i]*dimx*stages + nq_x[i]*dimx,
 		&blockidx[i], &blocksize[i], work+stid*lwork, &nrhs, &mrhs,
 		&nq_u[i], &nq_x[i], rkb, &hstep,
		&Q1[stid][blockidx[i]*dimx], (iQ2==0?0:&Q2[stid][(blockidx[i]-blocksize[0])*dimx]), &R[stid][blockidx[i]*3], &err, tempD, &iQ2, mergecholscheduler[i]);

	delete[] tempD;

	if(err!=0)
	{
		return CHOLFAILED;
	}

	return SUCCESS;

}

Kernel::T_ERR Kernel::mergeBlocksModule(unsigned int stid)
{

	int err = 0;

	int mrhs = Bcomputed[stid]==1?3:1;
	int nrhs = xsize+hsize;

	for (int i=0; i<N; i++)
	{
		int iQ2=(i>0?1:0);
		int cnidx = dimx*i;

		evdn_(&dimu, &dimx, &stages,
     		&blockidx[i], &blocksize[i], work+stid*lwork, &nrhs, &mrhs,
			&nq_u[i], &nq_x[i], rkb, &hstep,
			&Q1[stid][blockidx[i]*dimx], (iQ2==0?0:&Q2[stid][(blockidx[i]-blocksize[0])*dimx]), &R[stid][blockidx[i]*3], &iQ2,
			&DN[stid][dimx*dimx*i], (iQ2==0?0:&DN[stid][dimx*dimx*(i-1)]), (iQ2==0?0:&DN[stid][N*dimx*dimx + dimx*dimx*(i-1)]), &blockidx[N+1], &cnidx);
	}

	if (nq_x[N] > 0)
	{
		/**********************
		 * evbxrl			  *
		 **********************/
		int isq = nr_x[N]>0?1:0;
		int cnidx = blockidx[N+1]+dimx*(N-1);

		//TODO: workarray...
		double * temp = new double[dimx*nq_x[N]];
		double * D = new double [nq_x[N]*nq_x[N]];
		evbxrl_(&dimx,  sigx+stid*size_sigx+idx_sigx[N]-1, &isq, &nrhs,
     		&mrhs, work+stid*lwork, dqx+stid*size_dqx+idx_dq_x[N], &nq_x[N], &blockidx[N],
     		&Q1[stid][blockidx[N]*dimx], &R[stid][blockidx[N]*3],
     		temp, D,&DN[stid][dimx*dimx*(N-1)], &cnidx, &err, sigxcholscheduler[N]);
		delete[] temp;
		delete[] D;
	} else
	{
		/**********************
		 * evbxl			  *
		 **********************/
		int isq = nr_x[N]>0?1:0;
		int cnidx = blockidx[N+1]+dimx*(N-1);

		evbxl_(&dimx, sigx+stid*size_sigx+idx_sigx[N]-1, &isq, &nrhs, &mrhs,
			work+stid*lwork, &blockidx[N], &Q1[stid][blockidx[N]*dimx],
			&R[stid][blockidx[N]*3],&DN[stid][dimx*dimx*(N-1)], &cnidx, &err, sigxcholscheduler[N]);
	}
	if(err!=0)
	{
		return CHOLFAILED;
	}

	/**********************
	 * evdncn			  *
	 **********************/
	chlsch_(&N, &dimx, &mrhs, &nrhs, &blockidx[N+1], DN[stid],
	work+stid*lwork, &err, merge_ns3, merge_ns6, merge_nbs, merge_sdiag,
	merge_ss3, merge_ss6, merge_sbs);
	if(err!=0)
	{
		return CHOLFAILED;
	}

	for (int i=0; i<N+1; i++)
	{
		/**********************
		 * evsol			  *
		 **********************/
		int bidx = blockidx[i];
		int ynidx = blockidx[N+1]+(i<N?i:i-1)*dimx;
		int iQ2=(i>0 && i<N)?1:0;

		evsol_(&N, &dimx, &bidx, &ynidx, &blocksize[i], work+stid*lwork,
			&nrhs, &mrhs, &Q1[stid][blockidx[i]*dimx],
			(iQ2==0?0:&Q2[stid][(blockidx[i]-blocksize[0])*dimx]),
			 &R[stid][blockidx[i]*3], &iQ2, P, lgsb[stid]);
	}

	return SUCCESS;	
}

Kernel::T_ERR Kernel::coreModule(unsigned int stid)
{
	int err = 0;
	int nrhs = xsize+hsize;

	/**********************
	 * evcx				  *
	 **********************/
	if (Bcomputed[stid])
	{
		evcx_(&xsize, &hsize, y_b[stid], s_b[stid], xi+stid, lgsb[stid], cx1+stid, cx2+stid, &err);
		if(err!=0)
		{
			return GAUSSFAILED;
		}
	}

	/**********************
	 * evpxszk				  *
	 **********************/
	for (int i=0; i<N; i++)
	{
		evpxsz_(&dimu, &g[stid][idx_ub_u[i]],
	     &nlb_u[i], &nub_u[i], &(plb_u[val_idx_lb_u[i]-1]), &(pub_u[val_idx_ub_u[i]-1]),
	     &nr_u[i], &dru[stid][idx_dr_u[i]],
	     &lgsb[stid][i*dimu], &p_s[stid][idx_ub_u[i]], &p_z[stid][idx_ub_u[i]],
	     &vec_s[stid][idx_ub_u[i]], &vec_z[stid][idx_ub_u[i]], mu+stid,
	     cx1+stid, cx2+stid, &lgsb[stid][nrhs+i*dimu], &lgsb[stid][2*(nrhs)+i*dimu], Bcomputed+stid);
	}
	for (int i=0; i<N+1; i++)
	{
		evpxsz_(&dimx, &g[stid][idx_ub_x[i]],
	     &nlb_x[i], &nub_x[i], &(plb_x[val_idx_lb_x[i]-1]), &(pub_x[val_idx_ub_x[i]-1]),
	     &nr_x[i], &drx[stid][idx_dr_x[i]],
	     &lgsb[stid][N*dimu + i*dimx], &p_s[stid][idx_ub_x[i]], &p_z[stid][idx_ub_x[i]],
	     &vec_s[stid][idx_ub_x[i]], &vec_z[stid][idx_ub_x[i]], mu+stid,
	     cx1+stid, cx2+stid, &lgsb[stid][nrhs+N*dimu + i*dimx], &lgsb[stid][2*(nrhs)+N*dimu + i*dimx], Bcomputed+stid);
	}
	evpy_(&hsize, &lgsb[stid][xsize], cx1+stid, cx2+stid,  &lgsb[stid][(nrhs)+xsize],
		&lgsb[stid][2*(nrhs)+xsize], Bcomputed+stid);
	if (Bcomputed[stid])
	{
		for (int i=0; i<N; i++)
		{
			int idx = N*dimu + (N+1)*dimx+ i*dimk;
			evpxszk_(&dimk, &lgsb[stid][idx], cx1+stid, cx2+stid, 
			 &lgsb[stid][nrhs+idx], &lgsb[stid][2*(nrhs)+idx]);

		}				
	}

	/**********************
	 * evmer0			  *
	 **********************/
	phi0[stid]=0.0;
	dphi0[stid]=0.0;
	for (int i=0; i<N; i++)
	{
		evmer_(&dimu, &g[stid][idx_ub_u[i]], &h[stid][idx_q_u[i]],
			&nlb_u[i], &nub_u[i], &(plb_u[val_idx_lb_u[i]-1]), &(pub_u[val_idx_ub_u[i]-1]),
			&nr_u[i], &dru[stid][idx_dr_u[i]], &nq_u[i], dqu+stid*size_dqu+idx_dq_u[i],
			&lgsb[stid][i*dimu], &p_s[stid][idx_ub_u[i]],  &vec_s[stid][idx_ub_u[i]],
			mu+stid, phi0+stid, dphi0+stid);
	}
	for (int i=0; i<N+1; i++)
	{
		evmer_(&dimx, &g[stid][idx_ub_x[i]], &h[stid][idx_q_x[i]],
			&nlb_x[i], &nub_x[i], &(plb_x[val_idx_lb_x[i]-1]), &(pub_x[val_idx_ub_x[i]-1]),
			&nr_x[i], &drx[stid][idx_dr_x[i]], &nq_x[i], dqx+stid*size_dqx+idx_dq_x[i],
			&lgsb[stid][N*dimu + i*dimx], &p_s[stid][idx_ub_x[i]],  &vec_s[stid][idx_ub_x[i]],
			mu+stid, phi0+stid, dphi0+stid);
	}
	for (int i=0; i<N; i++)
	{
		evmeru_(&dimu, &dimx, &dimk, &stages, &h[stid][idx_hk[i*stages]],
			&lgsb[stid][i*dimu], dfu+stid*size_dfu+i*dimk*dimu, dphi0+stid);
		evmerx_(&dimx, &dimk, &stages, &h[stid][idx_hk[i*stages]],
			&lgsb[stid][N*dimu + i*dimx], &lgsb[stid][N*dimu + (N+1)*dimx + i*dimk],
			dfx+stid*size_dfx+i*dimk*dimx, &dfxA[stid][i*dimk*dimk],
     		rkb, &hstep, phi0+stid, dphi0+stid);

	}			

	/**********************
	 * stepsz			  *
	 **********************/
	stepsz_(&xsize, &gsize, &hsize, alpha_s+stid, alpha_z+stid, dphi0+stid,
			lgsb[stid], p_s[stid], p_z[stid], &lgsb[stid][xsize], minustau+stid, vec_s[stid], vec_z[stid]);

	int again=1;
	while (again)
	{
		/**********************
		 * evmeralpha		  *
		 **********************/
		merstp_(&xsize,&gsize,vec_x[stid],lgsb[stid],vec_s[stid],p_s[stid],alpha_s+stid,vec_x_tmp[stid],vec_s_tmp[stid]);

		int lwork=dimx+dimx*dimx+dimx*dimu;
		int dograd=0;
		int err=0;
		//TODO: work array...
		double * work = new double [lwork];

		evrest_(vec_x_tmp[stid], par, &t0, &N, &hstep,
	     &dimx, &dimu, &xsize, &gsize, &hsize,
	     &stages, rkA, rkb, rkc,
	     func_rhs, func_rhsx, func_rhsu, 
	     func_r_u, func_dr_u, func_r_x, func_dr_x,
	     func_q_u, func_dq_u, func_q_x, func_dq_x,
	     nub_u, nlb_u, nub_x, nlb_x,
	     val_idx_ub_u, val_idx_lb_u, val_idx_ub_x, val_idx_lb_x,
	     pub_u, plb_u, pub_x, plb_x,
	     ub_u, lb_u, ub_x[stid], lb_x[stid],
	     nr_u, nr_x, nq_u, nq_x,
	     idx_dr_u, idx_dr_x, idx_dq_u, idx_dq_x,
	     idx_ub_u, idx_lb_u, idx_r_u, idx_ub_x, idx_lb_x, idx_r_x,
	     idx_q_u, idx_q_x,
	     idx_hk, idx_hx,
	     g[stid], h[stid], dfu+stid*size_dfu, dfx+stid*size_dfx, dfxA[stid], dqu+stid*size_dqu, dqx+stid*size_dqx, dru[stid], drx[stid],
	     work, &lwork, &dograd, &err);

		delete[] work;

		phia[stid]=0.0;
		for (int i=0; i<N; i++)
		{
			evme2_(&dimu, &g[stid][idx_ub_u[i]], &h[stid][idx_q_u[i]],
				&nlb_u[i], &nub_u[i], &(plb_u[val_idx_lb_u[i]-1]), &(pub_u[val_idx_ub_u[i]-1]),
				&nr_u[i], &nq_u[i],
				&vec_s_tmp[stid][idx_ub_u[i]], mu+stid, phia+stid);
		}
		for (int i=0; i<N+1; i++)
		{
			evme2_(&dimx, &g[stid][idx_ub_x[i]], &h[stid][idx_q_x[i]],
				&nlb_x[i], &nub_x[i], &(plb_x[val_idx_lb_x[i]-1]), &(pub_x[val_idx_ub_x[i]-1]),
				&nr_x[i], &nq_x[i], 
				&vec_s_tmp[stid][idx_ub_x[i]], mu+stid, phia+stid);
		}
		for (int i=0; i<N; i++)
		{
			evmex2_(&dimx, &dimk, &h[stid][idx_hk[i*stages]], phia+stid);

		}			

		/**********************
		 * lnsrch			  *
		 **********************/
		lnsrch_(phi0+stid, phia+stid, dphi0+stid, alpha_s+stid, &eps, &again);
		if (again == -1)
		{
			return SMALLSTEPSIZE;
		}
	}

	/**********************
	 * evstp1			  *
	 **********************/
	for (int i=0; i<N && err==0; i++)
	{
		evstu1_(&dimu, &dimx, &dimk,
	     dfu+stid*size_dfu+dimu*dimk*i, &nq_u[i], dqu+stid*size_dqu+idx_dq_u[i],
	     &nr_u[i], &dru[stid][idx_dr_u[i]], &nlb_u[i], &nub_u[i],
	     &vec_x[stid][dimu*i], &lgsb[stid][dimu*i], &vec_s[stid][idx_ub_u[i]], &p_s[stid][idx_ub_u[i]],
	     &vec_z[stid][idx_ub_u[i]], &p_z[stid][idx_ub_u[i]],
	     &vec_y[stid][i*(dimk+dimx)], &lgsb[stid][xsize+i*(dimk+dimx)],
	     &vec_y[stid][idx_q_u[i]], &lgsb[stid][xsize+idx_q_u[i]], alpha_s+stid, alpha_z+stid,
	     &s_b_tmp[stid][dimu*i], &y_b_tmp[stid][dimu*i], &err);
	}
	for (int i=0; i<N+1 && err==0; i++)
	{
		int ti = i != N;
		evstx1_(&dimx, &dimk,
	     dfx+stid*size_dfx+dimx*dimk*i, &nq_x[i], dqx+stid*size_dqx+idx_dq_x[i],
	     &nr_x[i], &drx[stid][idx_dr_x[i]], &nlb_x[i], &nub_x[i],
	     &vec_x[stid][N*dimu+dimx*i], &lgsb[stid][dimu*N+dimx*i], &vec_s[stid][idx_ub_x[i]], &p_s[stid][idx_ub_x[i]],
	     &vec_z[stid][idx_ub_x[i]], &p_z[stid][idx_ub_x[i]],
	     &vec_y[stid][i*(dimk+dimx)],
	     &vec_y[stid][idx_q_x[i]], &lgsb[stid][xsize+idx_q_x[i]], alpha_s+stid, alpha_z+stid,
	     &s_b_tmp[stid][dimu*N+dimx*i], &y_b_tmp[stid][dimu*N+dimx*i], &ti, &err);
	}

	for (int i=0; i<N && err==0; i++)
	{
		evstk1_(&dimx, &dimk, &stages,
	     &dfxA[stid][dimk*dimk*i],
	     &vec_x[stid][N*dimu+dimx*(N+1)+dimk*i], &lgsb[stid][N*dimu+dimx*(N+1)+dimk*i], 
	     &vec_y[stid][i*(dimk+dimx)], alpha_s+stid,
	     &s_b_tmp[stid][N*dimu+dimx*(N+1)+dimk*i], &y_b_tmp[stid][N*dimu+dimx*(N+1)+dimk*i], &err);
	}
	if (err)
	{
		return FAILURE;
	}

	/**********************
	 * evrest    		  *
	 **********************/
	 evrest(stid, 1);

	/**********************
	 * evstp2    		  *
	 **********************/
	double sum=0.0;

	for (int i=0; i<N && err==0; i++)
	{
		evstu2_(&dimu, &dimx, &dimk,
	     dfu+stid*size_dfu+dimu*dimk*i, &nq_u[i], dqu+stid*size_dqu+idx_dq_u[i],
	     &nr_u[i], &dru[stid][idx_dr_u[i]], &nlb_u[i], &nub_u[i],
	     &vec_z[stid][idx_ub_u[i]], &vec_y[stid][i*(dimk+dimx)], &vec_y[stid][idx_q_u[i]],
	     &s_b_tmp[stid][dimu*i], &y_b_tmp[stid][dimu*i], &sum);

	}
	for (int i=0; i<N+1 && err==0; i++)
	{
		int ti = i != N;
		evstx2_(&dimx, &dimk,
	     dfx+stid*size_dfx+dimx*dimk*i, &nq_x[i], dqx+stid*size_dqx+idx_dq_x[i],
	     &nr_x[i], &drx[stid][idx_dr_x[i]], &nlb_x[i], &nub_x[i], 
	     &vec_z[stid][idx_ub_x[i]], &vec_y[stid][i*(dimk+dimx)], &vec_y[stid][idx_q_x[i]],
	     &s_b_tmp[stid][dimu*N+dimx*i], &y_b_tmp[stid][dimu*N+dimx*i], &ti, &sum);

	}

	for (int i=0; i<N && err==0; i++)
	{
		evstk2_(&dimx, &dimk, &stages,
	     &dfxA[stid][dimk*dimk*i],
	     &vec_y[stid][i*(dimk+dimx)],
	     &s_b_tmp[stid][N*dimu+dimx*(N+1)+dimk*i], &y_b_tmp[stid][N*dimu+dimx*(N+1)+dimk*i], &sum);
	}

	if (abs(sum)>1e-10 && false) /* TODO: BCompute wieder einschalten */
	{
		memcpy(s_b[stid], s_b_tmp[stid], sizeof(double)*xsize);
		memcpy(y_b[stid], y_b_tmp[stid], sizeof(double)*xsize);
		Bcomputed[stid]=1;
	} else
	{
		Bcomputed[stid]=0;
	}
	
	return SUCCESS;

}

Kernel::T_ERR Kernel::run(InitialGuess * init, double * xres)
{
	const int stid = 0;

	initModule(init, stid);
	bool oloop = true;
	while (oloop)
	{
		bool iloop = true;
		while (iloop)
		{
			T_ERR err;

			prepareSigRhsModule(stid);

			for (int i=0; i<N; i++)
			{
				err = innerBlockUModule(stid, i);
				if (err != SUCCESS) return err;

				err = innerBlockXModule(stid, i);
				if (err != SUCCESS) return err;
				
				err = innerBlockMergeModule(stid, i);
				if (err != SUCCESS) return err;
			}

			err = mergeBlocksModule(stid);
			if (err != SUCCESS) return err;

			err = coreModule(stid);
			if (err != SUCCESS) return err;

			iloop = innerLoopEndModule(stid, err);
			if (err != SUCCESS) return err;
		}
		oloop = outerLoopEndModule(stid);
	}

	if (steps[stid] >= maxsteps)
	{
		return MAXSTEPS;
	}

 	memcpy(xres, vec_x[stid], sizeof(double)*xsize);

	return SUCCESS;

}

void Kernel::getWarmstartData(double * z, double * y)
{
	memcpy(z, vec_z, gsize*sizeof(double));
	memcpy(y, vec_y, hsize*sizeof(double));
}

bool Kernel::getBFGSData(double * s_b, double * y_b)
{
	if (Bcomputed)
	{
		memcpy(s_b, this->s_b, xsize*sizeof(double));
		memcpy(y_b, this->y_b, xsize*sizeof(double));
	}
	return Bcomputed;
}

}