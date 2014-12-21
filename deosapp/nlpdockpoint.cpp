#include "nlpdockpoint.h"
#include "system.h"
#include <cstring>

using namespace Ipopt;

namespace DEOSApp
{

NLPDockpoint::NLPDockpoint(double * x, int N, int nqx, int nrx, double * par, double * targetN, double * vtargetN)
{
    this->N=N;
    this->nqx=nqx;
    this->nrx=nrx;
    this->par=par;
    this->x=x;
    this->vtargetN = vtargetN;
    this->targetN = targetN;
}

NLPDockpoint::~NLPDockpoint()
{
}

bool NLPDockpoint::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    n = 13;
    m = nqx+1+nrx+1;

    nnz_jac_g = (nqx+nrx)*13 + 4 + 6;

    nnz_h_lag = 2;

    index_style = FORTRAN_STYLE;

    return true;
}

bool NLPDockpoint::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
    for (int i=0; i<13; i++)
    {
        x_l[i] = -1e19;
        x_u[i] = 1e19;
    }

    // we have one equality constraint, so we set the bounds on this constraint
    // to be equal (and zero).
    for (int i=0; i<nqx+1; i++)
    {
        g_l[i] = g_u[i] = 0.0;
    }
    for (int i=nqx+1; i<nqx+1+nrx; i++)
    {
        g_l[i] = 0.0;
        g_u[i] = 1E19;
    }
    g_l[nqx+1+nrx] = 0.0;
    g_u[nqx+1+nrx] = 1E19;

    return true;
}

bool NLPDockpoint::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
    for (int i=0; i<13; i++)
    {
        x[i] = 0.0;
    }
    x[12] = 1.0;

    return true;
}

bool NLPDockpoint::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    double d=0.0;
    for (int i=0; i<12; i++)
    {
        d += x[i]*x[i];
    }
    d += (x[12]-1.0)*(x[12]-1.0);
    obj_value = d;
    return true;
}

bool NLPDockpoint::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    for (int i=0; i<12; i++)
    {
        grad_f[i] = 2.0*x[i];
    }
    grad_f[12] = 2.0*(x[12]-1.0);
    return true;
}

bool NLPDockpoint::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    int Np1=N+1;

    q_x_(&Np1, &N, (double*)x, par, g);
    g[nqx] = x[9]*x[9] + x[10]*x[10] + x[11]*x[11] + x[12]*x[12] - 1.0;
    r_x_(&Np1, &N, (double*)x, par, g+nqx+1);

    g[nqx+1+nrx] = (targetN[0] - x[0])*(x[3]-vtargetN[0]) + (targetN[1] - x[1])*(x[4]-vtargetN[1]) + (targetN[2] - x[2])*(x[5]-vtargetN[2]);
    return true;
}

bool NLPDockpoint::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
    if (values == NULL)
    {
        int cursor=0;

        for (int i=1; i<=13; i++)
        {   
            for (int j=1; j<=nqx; j++)
            {
                iRow[cursor] = j;
                jCol[cursor] = i;
                cursor++;
            }
        }

        for (int i=10; i<=13; i++)
        {
            iRow[cursor] = nqx+1;
            jCol[cursor] = i;
            cursor++;
        }

        for (int i=1; i<=13; i++)
        {   
            for (int j=1; j<=nrx; j++)
            {
                iRow[cursor] = nqx+1+j;
                jCol[cursor] = i;
                cursor++;
            }
        }
        for (int i=1; i<=6; i++)
        {
            iRow[cursor] = nqx+1+nrx+1;
            jCol[cursor] = i;
            cursor++;
        }
    } else
    {
        int Np1=N+1;

        dq_x_(&Np1, &N, (double*)x, par, values);
        values += nqx*13;

        values[0] = 2.0*x[9];
        values[1] = 2.0*x[10];
        values[2] = 2.0*x[11];
        values[3] = 2.0*x[12];
        values += 4;

        dr_x_(&Np1, &N, (double*)x, par, values);
        values += nrx*13;

        values[0] = -(x[3]-vtargetN[0]);
        values[1] = -(x[4]-vtargetN[1]);
        values[2] = -(x[5]-vtargetN[2]);
        values[3] = -(x[0]-targetN[0]);
        values[4] = -(x[1]-targetN[1]);
        values[5] = -(x[2]-targetN[2]);
    } //(targetN[0] - x[0])*(x[3]-vtargetN[0]) + (targetN[1] - x[1])*(x[4]-vtargetN[1]) + (targetN[2] - x[2])*(x[5]-vtargetN[2])

    return true;
}

bool NLPDockpoint::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
    return false;
}

void NLPDockpoint::finalize_solution(SolverReturn status,
                            Index n, const Number* x, const Number* z_L, const Number* z_U,
                            Index m, const Number* g, const Number* lambda,
                            Number obj_value,
                            const IpoptData* ip_data,
                            IpoptCalculatedQuantities* ip_cq)
{
    memcpy(this->x, x, sizeof(double)*13);
}

}

