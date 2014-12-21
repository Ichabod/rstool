#ifndef _DEOSAPP_DGLSYSTEM
#define _DEOSAPP_DGLSYSTEM

extern "C"
{
void rhs_(double * t, double * x, double * u, double * p, double * res);
void rhsx_(double * t, double * x, double * u, double * p, double * res);
void rhsu_(double * t, double * x, double * u, double * p, double * res);

void q_u_(int * i, int * N, double * u, double * p, double * res);
void dq_u_(int * i, int * N, double * u, double * p, double * res);
void nq_u_(int * N, double * p, int * res);
void r_u_(int * i, int * N, double * u, double * p, double * res);
void dr_u_(int * i, int * N, double * u, double * p, double * res);
void nr_u_(int * N, double * p, int * res);

void q_x_(int * i, int * N, double * x, double * p, double * res);
void dq_x_(int * i, int * N, double * x, double * p, double * res);
void nq_x_(int * N, double * p, int * res);
void r_x_(int * i, int * N, double * x, double * p, double * res);
void dr_x_(int * i, int * N, double * x, double * p, double * res);
void nr_x_(int * N, double * p, int * res);

void nbox_u_(int * N, double * p, int * lres, int * ures);
void pbox_u_(int * i, int * N, double * p, int * lres, int * ures);
void box_u_(int * i, int * N, double * p, double * lres, double * ures);

void nbox_x_(int * N, double * p, int * lres, int * ures);
void pbox_x_(int * i, int * N, double * p, int * lres, int * ures);
void box_x_(int * i, int * N, double * p, double * lres, double * ures);

void readnx_(int * dimx);
void readnu_(int * dimu);
void rdrks_(int * stages);
void readrk_(double * rkA, double * rkb, double * rkc);

void rotate_(double * x, double *res);
	
}

#endif
