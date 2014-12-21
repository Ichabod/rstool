#include "nlpmpc.h"
#include "system.h"
#include "rkstep.h"
#include <cstring>
#include <iostream>

using namespace Ipopt;
using namespace std;

namespace DEOSApp
{

NLPMPC::NLPMPC(double hstep, int mpcN, int N, double * target, double c_safety, double * par)
{
    this->hstep=hstep;
    this->mpcN=mpcN;
    this->N = N;
    this->target = target;
    this->par=par;
    this->c_safety = c_safety;

    rdrks_(&stages);

    rkA = new double[stages*stages];
    rkb = new double[stages];
    rkc = new double[stages];

    readrk_(rkA, rkb, rkc);

    rkwork = new double[(stages*(stages*13 + 2) + 1)*13];
    rkiwork = new int[stages*13];
    k= new double[13*stages];
    memset(k, 0, sizeof(double)*13*stages);

    g1 = new double[mpcN + 13];
    g2 = new double[mpcN + 13];

}

NLPMPC::~NLPMPC()
{
    delete[] rkA;
    delete[] rkb;
    delete[] rkc;
    delete[] rkwork;
    delete[] rkiwork;
    delete[] k;
    delete[] g1;
    delete[] g2;
}

void NLPMPC::setResultPtr(double * x, double * u)
{
    this->x = x;
    this->u = u;
}

void NLPMPC::setLastX(int idx, double * lastx)
{
    this->lastx = lastx;
    this->idx = idx;
}


void NLPMPC::rkstepWrapper(double hstep, double * x, double * u, double * par)
{
    int dim=13;
    double t= 0.0;
    double eps=1e-8;
    int ierr=0;

    rkstep_(&dim, &hstep, &t, x, x, u, par, &stages, k, rkA, rkb, rkc, &eps, rkwork, rkiwork, &ierr);
    if (ierr != 0) cout << "Fehler" << endl;
}


bool NLPMPC::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    n = mpcN*6 + 13;
    m = mpcN + 13;

    nnz_jac_g = m*n;

    nnz_h_lag = 0;

    index_style = FORTRAN_STYLE;

    return true;
}

bool NLPMPC::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
    for (int i=0; i<mpcN; i++)
    {
        x_l[i*6  ] = -par[10];
        x_l[i*6+1] = -par[10];
        x_l[i*6+2] = -par[10];
        x_l[i*6+3] = -par[11];
        x_l[i*6+4] = -par[11];
        x_l[i*6+5] = -par[11];
        x_u[i*6  ] = par[10];
        x_u[i*6+1] = par[10];
        x_u[i*6+2] = par[10];
        x_u[i*6+3] = par[11];
        x_u[i*6+4] = par[11];
        x_u[i*6+5] = par[11];
    }
    for (int i=0; i<13; i++)
    {
        x_l[6*mpcN + i] = -1e19;
        x_u[6*mpcN + i] = 1e19;
    }

    for (int i=0; i<mpcN; i++)
    {
        g_l[i] = 0.0;
        g_u[i] = 1E19;
    }

    for (int i=0; i<13; i++)
    {
        g_l[mpcN + i] = 0.0;
        g_u[mpcN + i] = 0.0;
    }

    return true;
}

bool NLPMPC::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
    memset(x, 0, sizeof(double)*6*mpcN);
    memcpy(x + 6*mpcN, lastx, sizeof(double)*13);

    return true;
}

bool NLPMPC::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    double x_ode[13];
    memcpy(x_ode, x+6*mpcN, sizeof(double)*13);

    double d1=0.0;
    double d2=0.0;
    for (int i=0; i<mpcN; i++)
    {
        rkstepWrapper(hstep, x_ode, (double*)x, par);

        for (int j=3; j<9; j++)
        {
            d1 += x_ode[j]*x_ode[j];
        }
        for (int j=9; j<12; j++)
        {
            d2 += x_ode[j]*x_ode[j];
        }
        d2 += (x_ode[12]-1.0)*(x_ode[12]-1.0);
        x += 6;
    }
    obj_value = d1+50.0*d2;

    return true;
}

bool NLPMPC::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    Number f1;
    Number f2;

    double * tmpx = (double*)x;

    eval_f(n,tmpx,new_x,f1);
    for (int i=0; i<n; i++)
    {
        double d = tmpx[i];
        tmpx[i] += 1e-6;
        eval_f(n,tmpx,new_x,f2);
        tmpx[i] = d;

        grad_f[i] = (f2-f1)/1e-6;
    }

    return true;
}

bool NLPMPC::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    double x_ode[13];
    memcpy(x_ode, x+6*mpcN, sizeof(double)*13);
    double * tmptarget = target+(idx-mpcN)*3;
    double * tmpx = (double*)x;
 
    for (int i=0; i<mpcN; i++)
    {
        double dx = x_ode[0] - tmptarget[0];
        double dy = x_ode[1] - tmptarget[1];
        double dz = x_ode[2] - tmptarget[2];

        *g = dx*dx + dy*dy + dz*dz - c_safety*c_safety;
        g++;

        rkstepWrapper(hstep, x_ode, tmpx, par);
        tmptarget += 3;
        tmpx+= 6;
    }

    for (int i=0; i<13; i++)
    {
        g[i] = x_ode[i] - lastx[i];
    }

    return true;
}

bool NLPMPC::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
    if (values == NULL)
    {
        int cursor=0;

        for (int i=1; i<=n; i++)
        {   
            for (int j=1; j<=m; j++)
            {
                iRow[cursor] = j;
                jCol[cursor] = i;
                cursor++;
            }
        }
    } else
    {
        int cursor=0;
        double * tmpx = (double*)x;
        eval_g(n,tmpx, new_x, m, g1);
        for (int i=0; i<n; i++)
        {
            double d = tmpx[i];
            tmpx[i] += 1e-6;
            eval_g(n,tmpx, new_x, m, g2);
            tmpx[i] = d;

            for (int j=0; j<m; j++)
            {
                values[cursor++] = (g2[j]-g1[j])/1e-6;
            }
        }

    }

    return true;
}

bool NLPMPC::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
    return false;
}

void NLPMPC::finalize_solution(SolverReturn status,
                            Index n, const Number* x, const Number* z_L, const Number* z_U,
                            Index m, const Number* g, const Number* lambda,
                            Number obj_value,
                            const IpoptData* ip_data,
                            IpoptCalculatedQuantities* ip_cq)
{
    memcpy(this->x, x+mpcN*6, sizeof(double)*13);
    for (int i=0; i<mpcN-1; i++)
    {
         rkstepWrapper(hstep, this->x, (double*)(x+i*6), par);
    }
    memcpy(this->u, x+(mpcN-1)*6, sizeof(double)*6);
}

}

