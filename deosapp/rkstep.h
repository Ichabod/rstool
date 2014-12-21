#ifndef _DEOSAPP_RKSTEP
#define _DEOSAPP_RKSTEP

extern "C"
{

void rkstep_(int * dim, double * h, double * t, double * x, double * xout, double * u, double * p, int * s, double * k, double * rkA, double * rkb, double * rkc, double * eps, double * work, int * iwork, int * ierr);
void newslv_(int * dim, double * h, double * t, double * x, double * u, double * p, double * k, int * s, double * eps, double * rkA, double * rkc, double * df, double * fres, double * fres2, int * ip, double * rktmp, int * ier);
void evrkf_(int * dim, double * h, double * t, double * x, double * u, double * p, double * k, int * s, double * rkA, double * rkc, double * res, double * tmp);
void evrkfu_(int * dim, int * dimu, double * h, double * t, double * x, double * u, double * p, double * k, int * s, double * rkA, double * rkc, double * res, double * tmp, double * tmp2);
void evrkfx_(int * dim, double * h, double * t, double * x, double * u, double * p, double * k, int * s, double * rkA, double * rkc, double * res, double * tmp, double * tmp2);

void sol_(int * n, int * ndim, double * A, double * b, int * ip);

}

#endif
