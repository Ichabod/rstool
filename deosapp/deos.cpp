#include "deos.h"
#include "system.h"
#include "rkstep.h"
#include "../ode/dopri5.h"
#include "../ode/exception.h"
#include "nlpdockpoint.h"
#include "nlpmpc.h"
#include "nlpmpcsparse.h"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>

using namespace std;
using namespace Ipopt;

namespace DEOSApp
{

DEOS::DEOS(int N, double hstep, double c_rx, double c_G, double c_Ms,
            double c_js1, double c_js2, double c_js3,
            double c_dsx, double c_dsy, double c_dsz, double c_vmax, double c_mmax, double c_safety)
{
    this->N = N;
    this->hstep = hstep;
    
    this->c_G=c_G;
    this->c_rx = c_rx;

    par=new double[12 + (N+1)*3 + (N+1)*3 + (N+1)*3];
    par[ 0] = sqrt(c_G/pow(c_rx,3.0));
    par[ 1] = 0.0;
    par[ 2] = this->c_Ms = c_Ms;
    par[ 3] = this->c_js1 = c_js1;
    par[ 4] = this->c_js2 = c_js2;
    par[ 5] = this->c_js3 = c_js3;
    par[ 6] = this->c_dsx = c_dsx;
    par[ 7] = this->c_dsy = c_dsy;
    par[ 8] = this->c_dsz = c_dsz;
    par[ 9] = this->c_safety = c_safety;
    par[10] = this->c_vmax = c_vmax;
    par[11] = this->c_mmax = c_mmax;

    vtargetN = new double[3];
    target = par+12;
    Rtd = target+(N+1)*3;
    dRtd = Rtd+(N+1)*3;

    ipoptapp = IpoptApplicationFactory();//new IpoptApplication();
    ipoptapp->Initialize();

    ipoptapp->Options()->SetStringValue ( "sb", "yes" );
    ipoptapp->Options()->SetIntegerValue ( "print_level", 0 );
    ipoptapp->Options()->SetIntegerValue ( "limited_memory_max_history", 10 );
    ipoptapp->Options()->SetStringValue ( "limited_memory_update_type", "bfgs" );
    ipoptapp->Options()->SetStringValue ( "limited_memory_initialization", "scalar1" );
    ipoptapp->Options()->SetNumericValue ( "limited_memory_init_val", 1 );
    ipoptapp->Options()->SetNumericValue ( "limited_memory_init_val_max", 1E8 );
    ipoptapp->Options()->SetNumericValue ( "limited_memory_init_val_min", 1E-8 );
    ipoptapp->Options()->SetIntegerValue ( "limited_memory_max_skipping", 5 );
    ipoptapp->Options()->SetStringValue ( "hessian_approximation", "limited-memory" );
    ipoptapp->Options()->SetNumericValue ( "tol", 1E-6 );
//    ipoptapp->Options()->SetStringValue ( "linear_solver", "mumps" );
   

    nqx = new int[N+1];
    nrx = new int[N+1];

    nq_x_(&N,par,nqx);
    nr_x_(&N,par,nrx);

    if (nqx[N] > 0)
    {
        dqxNtemp = new double[nqx[N]*13];   
    } else
    {
        dqxNtemp = 0;
    }

    if (nrx[N] > 0)
    {
        drxNtemp = new double[nrx[N]*13];   
    } else
    {
        drxNtemp = 0;
    }

    mpcN=max(min(N/4, 10),5);;
    mpcstart = new double[13];

    rdrks_(&stages);

    rkA = new double[stages*stages];
    rkb = new double[stages];
    rkc = new double[stages];

    readrk_(rkA, rkb, rkc);

    rkwork = new double[(stages*(stages*13 + 2) + 1)*13];
    rkiwork = new int[stages*13];
    k= new double[13*stages];
    memset(k, 0, sizeof(double)*13*stages);

    mpcsolution = new double[mpcN*6 + mpcN*13];
}

DEOS::~DEOS()
{
    delete[] par;
    delete[] nqx;
    delete[] nrx;
    if (dqxNtemp != 0) delete[] dqxNtemp;
    if (drxNtemp != 0) delete[] drxNtemp;
    delete[] mpcstart;
    delete[] rkA;
    delete[] rkb;
    delete[] rkc;
    delete[] rkwork;
    delete[] rkiwork;
    delete[] k;
    delete[] vtargetN;
    delete[] mpcsolution;

}

void DEOS::ODEwrapper( int *, double * t, double * y, double * dy, double * rpar, int * )
{
    rhs_(t, y, rpar+6, rpar, dy);
}


bool DEOS::generateTarget(double * xt, double c_Mt, double c_jt1, double c_jt2, double c_jt3, double c_dtx, double c_dty, double c_dtz)
{
    double tpar[] = { sqrt(c_G/pow(c_rx,3.0)), 0.0, c_Mt, c_jt1, c_jt2, c_jt3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double dt[] = {c_dtx, c_dty, c_dtz};
    double Rmat[9];
    double Rtw[3];

    double x[13];
    memcpy(x, xt, sizeof(double)*13);

    try
    {
        ODE::DoPri5 ode(ODEwrapper, 13);
        ode.setTol(1e-8, 1e-8);
        ode.init(0.0, x);

        double * Rtditem = Rtd;
        double * dRtditem = dRtd;
        for (int i=0; i<N+1; i++)
        {
            if (i > 0)
            {
                ode.calc((double)i*hstep, tpar);
            }
            memcpy(target+(i*3), x, sizeof(double)*3);
            if (i==N) memcpy(vtargetN, x+3, sizeof(double)*3);
            rotate_(x, Rmat);

            for (int j=0; j<3; j++)
            {
                double d1=0.0;
                double d2=0.0;
                for (int k=0; k<3; k++)
                {
                    double r = Rmat[k*3+j];
                    d1 += r*dt[k];
                    d2 += r*x[6+k];
                }
                Rtditem[j] = d1;
                Rtw[j] = d2;
            }
            for (int j=0; j<3; j++)
            {
                dRtditem[j] = Rtw[(j+1)%3]*Rtditem[(j+2)%3] - Rtw[(j+2)%3]*Rtditem[(j+1)%3];
            }
            Rtditem +=3;
            dRtditem +=3;
        }
        return true;
    } catch (ODE::Exception & e)
    {
        return false;
    }
}



bool DEOS::rkstepWrapper(double hstep, double * x, double * x_to, double * k, double * par)
{
    int dim=13;
    double t= 0.0;
    double eps=1e-8;
    int ierr=0;

    rkstep_(&dim, &hstep, &t, x,x_to, par+6, par, &stages, k, rkA, rkb, rkc, &eps, rkwork, rkiwork, &ierr);    
    return ierr==0;
}

bool DEOS::rkstepWrapper(double hstep, double * x, double * par)
{
    int dim=13;
    double t= 0.0;
    double eps=1e-8;
    int ierr=0;

    rkstep_(&dim, &hstep, &t, x,x, par+6, par, &stages, k, rkA, rkb, rkc, &eps, rkwork, rkiwork, &ierr);
    return ierr==0;
}

void DEOS::printVec(int step, double * x, int n)
{
    cout << "i = " << setw(2) << step << ": ";
    for (int i=0; i<n; i++)
    {
        cout << setw(9) << setprecision(5) << fixed << x[i];
    }
    cout << endl;
}

int DEOS::generateInitialGuess(Algorithm::InitialGuess * & initialguessobj, double * poshint)
{
    ApplicationReturnStatus status;

    initialguessobj = new Algorithm::InitialGuess();
    double * initialguess = initialguessobj->x = new double[Algorithm::InitialGuess::xsize];

    memset(initialguess, 0, sizeof(double)*(N*6 + (N+1)*13 + N*stages*13));

    /* Dockingpunkt berechnen */
    SmartPtr<NLPDockpoint> dock = new NLPDockpoint(mpcstart, N, nqx[N], nrx[N], par, target+N*3, vtargetN);
//    ipoptapp->Options()->SetStringValue ( "derivative_test", "first-order");
    ipoptapp->Options()->SetNumericValue ( "tol", 1E-9 );
    status = ipoptapp->OptimizeTNLP(dock);
    if (status < 0)
    {
        delete initialguessobj;
        initialguessobj=0;
        return -1;
    }
//    printVec(N, mpcstart, 13);
    memcpy(initialguess + N*6 + N*13, mpcstart, sizeof(double)*13);

    /* Startschätzung für MPC */
    memset(mpcsolution, 0, sizeof(double)*mpcN*6);
    double xtemp[13];
    memcpy(xtemp, mpcstart, sizeof(double)*13);
    double tpar[] = { sqrt(c_G/pow(c_rx,3.0)), 0.0, c_Ms, c_js1, c_js2, c_js3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i=0; i<mpcN; i++)
    {
        if (!rkstepWrapper(-hstep,xtemp,tpar))
        {
            delete initialguessobj;
            initialguessobj=0;
            return -2;
        }
        memcpy(mpcsolution+6*mpcN + (mpcN-1-i)*13, xtemp, sizeof(double)*13);
    }

    /* MPC */
    SmartPtr<NLPMPCSparse> mpc = new NLPMPCSparse(hstep, mpcN, N, target, mpcstart, poshint, mpcsolution, c_safety, par);


//    ipoptapp->Options()->SetIntegerValue ( "print_level", 5 );
//    ipoptapp->Options()->SetStringValue ( "expect_infeasible_problem", "yes");
    ipoptapp->Options()->SetIntegerValue ( "limited_memory_max_history", 5*mpcN );
//    ipoptapp->Options()->SetStringValue ( "derivative_test", "first-order");
    ipoptapp->Options()->SetStringValue ("nlp_scaling_method", "user-scaling");
    ipoptapp->Options()->SetNumericValue ( "tol", 1E-6 );

    for (int i=N-1; i>=mpcN; i--)
    {
        if (i==N-1)
        {
            status = ipoptapp->OptimizeTNLP(mpc);
            ipoptapp->Options()->SetStringValue ("warm_start_init_point", "yes" );
        } else
        {
            status = ipoptapp->ReOptimizeTNLP(mpc);
        }
        if (status < 0)
        {
            delete initialguessobj;
            initialguessobj=0;
            return -3;
        }
        mpc->shift(initialguess+i*6, initialguess+N*6 + i*13, initialguess + N*6 + (N+1)*13 + i*stages*13, true);
 //       printVec(i, initialguess+N*6 + i*13, 13);
    }

    status = ipoptapp->ReOptimizeTNLP(mpc);
    if (status < 0)
    {
        delete initialguessobj;
        initialguessobj=0;
        return -3;
    }
    memcpy(initialguess, mpcsolution, sizeof(double)*mpcN*6);
    memcpy(initialguess+N*6, mpcsolution+mpcN*6, sizeof(double)*mpcN*13);
    for (int i=mpcN-1; i>=0; i--)
    {
        memcpy(tpar+6, mpcsolution+i*6, sizeof(double)*6);
        rkstepWrapper(hstep, mpcsolution+mpcN*6+i*13, mpcsolution+mpcN*6+i*13, initialguess + N*6 + (N+1)*13 + i*stages*13, tpar);
//        printVec(i, initialguess+N*6 + i*13, 13);
    }

    return 0;
}




}