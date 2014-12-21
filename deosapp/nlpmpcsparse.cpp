#include "nlpmpcsparse.h"
#include "system.h"
#include "rkstep.h"
#include <cstring>
#include <iostream>
#include <cmath>

using namespace Ipopt;
using namespace std;

#define CTRLWEIGHT 1.0e-4
#define ROTCTRLWEIGHT 1.0e-4
#define POSWEIGHT 0.0001
#define SPEEDWEIGHT 10.0
#define ROTSPEEDWEIGHT 1000.0
#define ROTWEIGHT 10000.0

namespace DEOSApp
{

NLPMPCSparse::NLPMPCSparse(double hstep, int mpcN, int N, double * target, double * lastx, double * poshint, double * solution, double c_safety, double * par)
{
    this->hstep=hstep;
    this->mpcN=mpcN;
    this->N = N;
    this->target = target;
    this->par=par;
    this->c_safety = c_safety;
    this->poshint = poshint;

    rdrks_(&stages);

    rkA = new double[stages*stages];
    rkb = new double[stages];
    rkc = new double[stages];

    readrk_(rkA, rkb, rkc);

    rkwork = new double[(stages*(stages*13 + 2) + 1)*13];
    rkiwork = new int[stages*13];
    workfres = rkwork + stages*13*stages*13;
    workfres2 = workfres + stages*13;
    workrktmp = workfres2 + stages*13;

    k= new double[mpcN*13*stages];
    df = new double[mpcN*13*stages*13*stages];
    ip = new int[mpcN*13*stages];
    dfu = new double[13*stages*6];
    dfx = new double[13*stages*13];

    this->lastx = new double[13];
    memcpy(this->lastx, lastx, sizeof(double)*13);
    idx=N;

    this->solution = solution;
    z_L = new double[mpcN*(6+13)];
    z_U = new double[mpcN*(6+13)];
    for (int i=0; i<mpcN*(6+13); i++)
    {
        z_L[i] = z_U[i] = 0.1;
    }
    lambda = new double[mpcN*(13+1)];
    for (int i=0; i<mpcN*(13+1); i++)
    {
        lambda[i] = 1.0;
    }

    memset(k, 0, sizeof(double)*mpcN*13*stages);

    obj_scaling=1.0/(100.0 * (double)mpcN);
    restart = false;

}

NLPMPCSparse::~NLPMPCSparse()
{
    delete[] rkA;
    delete[] rkb;
    delete[] rkc;
    delete[] rkwork;
    delete[] rkiwork;
    delete[] k;
    delete[] df;
    delete[] ip;
    delete[] dfu;
    delete[] dfx;
    delete[] lastx;
    delete[] z_L;
    delete[] z_U;
    delete[] lambda;
}

void NLPMPCSparse::recomputeRK(const double * xvec, const double * uvec)
{
    int dim=13;
    double t=0.0;
    const int kskip=stages*13;
    const int dfskip=kskip*kskip;

    double * tdf = df;
    double * tk = k;
    int * tip = ip;
    int ier = 0;
    double * tmpxvec = (double*)xvec;
    double * tmpuvec = (double*)uvec;

    for (int i=0; i<mpcN; i++)
    {
        double eps=1e-8;
        newslv_(&dim, &hstep, &t, tmpxvec, tmpuvec, par, tk, &stages, &eps, rkA, rkc, tdf, workfres, workfres2, tip, workrktmp, &ier);
        if (ier != 0) cout << "FEHLER: " << ier << endl;

        tmpxvec += 13;
        tmpuvec += 6;
        tk += kskip;
        tdf += dfskip;
        tip += kskip;
    }
}

void NLPMPCSparse::shift(double * u, double * x, double * k, bool domemshift)
{
    int dim=13;
    double t=0.0;
    double eps=1e-8;
    int ierr=0;
    double tempx[13];

    memcpy(x, solution+mpcN*6+(mpcN-1)*13, sizeof(double)*13);
    memcpy(u, solution+(mpcN-1)*6, sizeof(double)*6);

    rkstep_(&dim, &hstep, &t, x, tempx, u, par, &stages, k, rkA, rkb, rkc, &eps, rkwork, rkiwork, &ierr);
    if (ierr != 0) cout << "Fehler" << endl;

    memcpy(lastx, x, sizeof(double)*13);

    if (domemshift)
    {
        memmove(solution+6, solution, sizeof(double)*(mpcN-1)*6);
        memmove(solution+mpcN*6+13, solution+mpcN*6, sizeof(double)*(mpcN-1)*13);
        for (int i=0; i<mpcN*(6+13); i++) solution[i] += 1e-3;

        memmove(z_L+6, z_L, sizeof(double)*(mpcN-1)*6);
        memmove(z_U+6, z_U, sizeof(double)*(mpcN-1)*6);
        memmove(z_L+mpcN*6+13, z_L+mpcN*6, sizeof(double)*(mpcN-1)*13);
        memmove(z_U+mpcN*6+13, z_U+mpcN*6, sizeof(double)*(mpcN-1)*13);
        memmove(lambda+13, lambda, sizeof(double)*(mpcN-1)*13);
        memmove(lambda+mpcN*13+1, lambda+mpcN*13, sizeof(double)*(mpcN-1));
    }

    idx--;
}


void NLPMPCSparse::rkstepWrapper(double hstep, double * x, double * u, double * par)
{
    int dim=13;
    double t= 0.0;
    double eps=1e-8;
    int ierr=0;

    rkstep_(&dim, &hstep, &t, x, x, u, par, &stages, k, rkA, rkb, rkc, &eps, rkwork, rkiwork, &ierr);
    if (ierr != 0) cout << "Fehler" << endl;
}


bool NLPMPCSparse::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    n = mpcN*(6 + 13);
    m = mpcN*(1 + 13);

    nnz_jac_g = mpcN*( 13*6 + 13*13 + 3) + (mpcN-1)*13;

    nnz_h_lag = 0;

    index_style = FORTRAN_STYLE;

    return true;
}

bool NLPMPCSparse::get_bounds_info(Index n, Number* x_l, Number* x_u,
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
    for (int i=0; i<13*mpcN; i++)
    {
        x_l[6*mpcN + i] = -1e19;
        x_u[6*mpcN + i] = 1e19;
    }

    for (int i=0; i<13*mpcN; i++)
    {
        g_l[i] = 0.0;
        g_u[i] = 0.0;
    }

    for (int i=0; i<mpcN; i++)
    {
        g_l[13*mpcN + i] = 0.0;
        g_u[13*mpcN + i] = 1E19;
    }

    return true;
}

bool NLPMPCSparse::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
    
    memcpy(x, solution, sizeof(double)*mpcN*(6+13));
    if (init_z)
    {
        memcpy(z_L, this->z_L, sizeof(double)*n);
        memcpy(z_U, this->z_U, sizeof(double)*n);
    }
    if (init_lambda)
    {
        memcpy(lambda, this->lambda, sizeof(double)*m);
    }

    return true;
}

bool NLPMPCSparse::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    new_x = true;
    if (new_x) recomputeRK(x+mpcN*6, x);

    double * x_ode = (double*)x+6*mpcN;
    double * u_ode = (double*)x;

    double d=0.0;
    for (int i=0; i<mpcN; i++)
    {
        for (int j=0; j<3; j++)
        {
            d += CTRLWEIGHT*u_ode[j]*u_ode[j];
        }
        for (int j=3; j<6; j++)
        {
            d += ROTCTRLWEIGHT*u_ode[j]*u_ode[j];
        }
        for (int j=0; j<3; j++)
        {
            d += POSWEIGHT*(x_ode[j]-poshint[j])*(x_ode[j]-poshint[j]);
        }
        for (int j=3; j<6; j++)
        {
            d += SPEEDWEIGHT*x_ode[j]*x_ode[j];
        }
        for (int j=6; j<9; j++)
        {
            d += ROTSPEEDWEIGHT*x_ode[j]*x_ode[j];
        }
        for (int j=9; j<12; j++)
        {
            d += ROTWEIGHT*x_ode[j]*x_ode[j];
        }
        d +=  ROTWEIGHT*(x_ode[12]-1.0)*(x_ode[12]-1.0);
        x_ode += 13;
        u_ode += 6;
    }
    obj_value = d;

    //cout << "f: " << obj_value << endl;

    return true;
}

bool NLPMPCSparse::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    new_x = true;
    if (new_x) recomputeRK(x+mpcN*6, x);

    double * x_ode = (double*)x+6*mpcN;
    double * u_ode = (double*)x;
    double * grad_u = grad_f;
    double * grad_x = grad_f+6*mpcN;

    for (int i=0; i<mpcN; i++)
    {
        for (int j=0; j<3; j++)
        {
            grad_u[j] = 2.0*CTRLWEIGHT*u_ode[j];
        }
        for (int j=3; j<6; j++)
        {
            grad_u[j] = 2.0*ROTCTRLWEIGHT*u_ode[j];
        }
        for (int j=0; j<3; j++)
        {
            grad_x[j] = 2.0*POSWEIGHT*(x_ode[j]-poshint[j]);
        }
        for (int j=3; j<6; j++)
        {
            grad_x[j] = 2.0*SPEEDWEIGHT*x_ode[j];
        }
        for (int j=6; j<9; j++)
        {
            grad_x[j] = 2.0*ROTSPEEDWEIGHT*x_ode[j];
        }
        for (int j=9; j<12; j++)
        {
            grad_x[j] = 2.0*ROTWEIGHT*x_ode[j];
        }
        grad_x[12] = 2.0*ROTWEIGHT*(x_ode[12]-1.0);

        x_ode += 13;
        u_ode += 6;
        grad_x += 13;
        grad_u += 6;
    }

    return true;
}

bool NLPMPCSparse::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    if (new_x) recomputeRK(x+mpcN*6, x);

    double * x_ode = (double*)x+mpcN*6;
    double * g_rk = g;
    double * g_dist = g+mpcN*13;
    double * tmptarget = target+(idx-mpcN)*3;
    double * tmpk = k;

    for (int i=0; i<mpcN; i++)
    {
        double * x2 = (i==mpcN-1 ? lastx : x_ode+13);

        double testx[13];
        memcpy(testx,x_ode,sizeof(double)*13);
        rkstepWrapper(hstep, testx, (double*)x+i*6, par);

        for (int j=0; j<13; j++)
        {
            double d = 0.0;
            double * rkk = tmpk;
            for (int j1=0; j1<stages; j1++)
            {
                d += rkb[j1]*rkk[j];
                rkk += 13;
            }
            g_rk[j] = x_ode[j] + hstep*d - x2[j];
        }
        double dx = x_ode[0] - tmptarget[0];
        double dy = x_ode[1] - tmptarget[1];
        double dz = x_ode[2] - tmptarget[2];
        g_dist[0] = dx*dx + dy*dy + dz*dz - c_safety*c_safety;

        g_rk += 13;
        g_dist += 1;
        tmptarget += 3;
        tmpk += 13*stages;
        x_ode += 13;
    }
 
    return true;
}

bool NLPMPCSparse::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{

    if (values == NULL)
    {
        int cursor=0;

        for (int i=0; i<mpcN; i++)
        {
            //du
            for (int j1=1; j1<=6; j1++)
            {
                for (int i1=1; i1<=13; i1++)
                {
                    iRow[cursor] = i*13 + i1;
                    jCol[cursor] = i*6 + j1;
                    cursor++;
                }
            }
            //dx
            for (int j1=1; j1<=13; j1++)
            {
                for (int i1=1; i1<=13; i1++)
                {
                    iRow[cursor] = i*13 + i1;
                    jCol[cursor] = mpcN*6 + i*13 + j1;
                    cursor++;
                }
            }
            //dist
            for (int j1=1; j1<=3; j1++)
            {
                iRow[cursor] = mpcN*13 + (i+1);
                jCol[cursor] = mpcN*6 + i*13 + j1;
                cursor++;
            }
        }
        // -ident
        for (int i=1; i<=(mpcN-1)*13; i++)
        {
            iRow[cursor] = i;
            jCol[cursor] = mpcN*6 + 13 + i;
            cursor++;
        }

    } else
    {
        if (new_x) recomputeRK(x+mpcN*6, x);

        double * x_ode = (double*)x+mpcN*6;
        double * u_ode = (double*)x;
        double * tmptarget = target+(idx-mpcN)*3;
        double * tmpk = k;
        double * tdf = df;
        int * tip = ip;

        const int ipskip = stages*13; 
        const int dfskip = ipskip*ipskip;
        int dim=13;
        int tdim=13*stages;
        int dimu=6;
        double t=0.0;
        double tmpdfu[13*6];
        double tmpdfx[13*13];

        for (int i=0; i<mpcN; i++)
        {
            /* DU */
            double * tdfu = dfu;

            evrkfu_(&dim, &dimu, &hstep, &t, x_ode, u_ode, par, tmpk, &stages, rkA, rkc, tdfu, workrktmp, tmpdfu);
            for (int l=0; l<6; l++)
            {
                sol_(&tdim, &tdim, tdf, tdfu, tip);
                for (int j=0; j<13; j++)
                {
                    double d = 0.0;
                    for (int j1=0; j1<stages; j1++)
                    {
                        d += rkb[j1]*tdfu[j1*13 + j];
                    }
                    *values = hstep*d;

                    values++;
                }
                tdfu += stages*13;
            }

            /* DX */
            double * tdfx = dfx;

            evrkfx_(&dim, &hstep, &t, x_ode, u_ode, par, tmpk, &stages, rkA, rkc, tdfx, workrktmp, tmpdfx);
            for (int l=0; l<13; l++)
            {
                sol_(&tdim, &tdim, tdf, tdfx, tip);
                for (int j=0; j<13; j++)
                {
                    double d = 0.0;
                    for (int j1=0; j1<stages; j1++)
                    {
                        d += rkb[j1]*tdfx[j1*13 + j];
                    }
                    d *= hstep;
                    if (l==j) d += 1.0;
                    *values = d;
                    values++;
                }
                tdfx += stages*13;
            }

            /* DIST */
            for (int j=0; j<3; j++)
            {
                *values = 2.0*(x_ode[j] - tmptarget[j]);
                values++;
            }

            tmptarget += 3;
            x_ode += 13;
            u_ode += 6;
            tdf += dfskip;
            tip += ipskip;
            tmpk += ipskip;
        }

        for (int i=0; i<(mpcN-1)*13; i++)
        {
            values[i] = -1.0;
        }


    }

    return true;
}

bool NLPMPCSparse::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
    return false;
}

void NLPMPCSparse::finalize_solution(SolverReturn status,
                            Index n, const Number* x, const Number* z_L, const Number* z_U,
                            Index m, const Number* g, const Number* lambda,
                            Number obj_value,
                            const IpoptData* ip_data,
                            IpoptCalculatedQuantities* ip_cq)
{
    memcpy(solution, x, sizeof(double)*n);
    memcpy(this->z_L, z_L, sizeof(double)*n);
    memcpy(this->z_U, z_U, sizeof(double)*n);
    memcpy(this->lambda, lambda, sizeof(double)*m);

    restart = true;
    obj_scaling = 1.0/(100.0*abs(obj_value));
}

bool NLPMPCSparse::get_scaling_parameters(Number& obj_scaling,
                                    bool& use_x_scaling, Index n,
                                    Number* x_scaling,
                                    bool& use_g_scaling, Index m,
                                    Number* g_scaling)
{
    use_x_scaling=false;
    use_g_scaling=false;

    obj_scaling=this->obj_scaling;
    if (restart) obj_scaling *= 4*mpcN;

    return true;
}

}

