typedef void(*T_RHS)(double *, double *, double *, double *, double *);
typedef void(*T_REST)(int *, int *, double *, double *, double *);

extern "C"
{

void evrest_(double * x, double * p, double * t0, int * N, double * hstep,
     int * dimx, int * dimu, int * xsize, int * gsize, int * hsize,
     int * s, double * A, double * b, double * c,
     T_RHS rhs, T_RHS rhsx, T_RHS rhsu,
     T_REST r_u, T_REST dr_u, T_REST r_x, T_REST dr_x,
     T_REST q_u, T_REST dq_u, T_REST q_x, T_REST dq_x,
     int * nub_u, int * nlb_u, int * nub_x, int * nlb_x,
     int * iub_u, int * ilb_u, int * iub_x, int * ilb_x,
     int * pub_u, int * plb_u, int * pub_x, int * plb_x,
     double * ub_u, double * lb_u, double * ub_x, double * lb_x,
     int * nr_u, int * nr_x, int * nq_u, int * nq_x,
     int * idr_u, int * idr_x, int * idq_u, int * idq_x,
     int * igub_u, int * iglb_u, int * ir_u, int * igub_x, int * iglb_x, int * ir_x,
     int * iq_u, int * iq_x,
     int * ihk, int * ihx,
     double * g, double * h, double * dfu, double * dfx, double * dfxA, double * dqu, double * dqx, double * dru, double * drx,
     double * work, int * lwork, int * dograd, int * err);

void evsig_(double * xi, double * vec_s, double * vec_z, int * N, int * dimx, int * dimu, int * gsize,
     int * nub_u, int * nlb_u, int * nub_x, int * nlb_x,
     int * pub_u, int * plb_u, int * pub_x, int * plb_x,
     int * nr_u, int * nr_x,
     int * iub_u, int * ilb_u, int * ir_u, int * iub_x, int * ilb_x, int * ir_x,
     int * isigu, int * isigx, double * sigu, double * sigx,
     double * dru, double * drx);

void initip_(int * xsize, int * gsize, int * hsize, int * istart,
     double * mu, double * mtau, double * ceps, int * steps, int * osteps, int * cinf, 
     double * s, double * z, double * y);

void evrhs_(double * rhs, int * dimx, int * dimu, int * dimk, int * xsize, int * hsize, int * N,
     int * nub_u, int * nlb_u, int * nub_x, int * nlb_x,
     int * pub_u, int * plb_u, int * pub_x, int * plb_x,
     int * nr_u, int * nr_x, double * r_u, double * r_x,
     int * nq_u, int * nq_x, double * q_u, double * q_x,
     double * g, double * h, double * dfu, double * dfx, double * dfxA,
     int * rks, double * rkb, double * hstep,
     double * y, double * z, double * s, double * mu,
     double * yb, double * sb, int * hasB, double * xi);

void evbu_(int * dimu, int * dimx, int * stages, double * sigu, int * isq, double * dfu, int * nrhs, 
     int * mrhs, double * rhs, int * bidx, double * y13, int * err, int * pinfo);

void evbx_(int * dimx, int * stages, double * sigx, int * isq, double * dfx, int * nrhs, 
     int * mrhs, double * rhs, int * bidx, double * y579, int * err, int * pinfo);

void evbxl_(int * dimx, double * sigx, int * isq, int * nrhs,
     int * mrhs, double * rhs, int * bidx, double * Q, double * R, double * DN, int * cNidx, int * err, int * pinfo);

void evbur_(int * dimu, int * dimx, int * stages, double * sigu, int * isq, double * dfu, int * nrhs,
     int * mrhs, double * rhs, double * dqu, int * nqu, int * bidx, double * y13, double * y24, double * D, double * temp, int * err, int * pinfo);

void evbxr_(int * dimx, int * stages, double * sigx, int * isq, double * dfx, int * nrhs,
     int * mrhs, double * rhs, double * dqx, int * nqx, int * bidx, double * y579, double * y6810, double * D, double * temp, int * err, int * pinfo);

void evbxrl_(int * dimx, double * sigx, int * isq, int * nrhs,
     int * mrhs, double * rhs, double * dqx, int * nqx, int * bidx,
     double * Q, double * R, double * temp, double * D, double * DN, int * cNidx, int * err, int * pinfo);

void evmrpr_(int * dimu, int * dimx, int * stages, double * dfu, double * dfx, double * dfxA, double * xi,
     double * y1, double * y3, double * y5, double * y7, double * y9, int * bidx, int * bsz, double * rhs, int * nrhs, int * mrhs,
     int * nqu, int * nqx, double * rkb, double * hstep,
     double * Q1, double * Q2, double * R, int * err, double * D, int * iQ2);

void evmrge_(int * dimu, int * dimx, int * stages, double * dfxA, double * xi,
     double * y1, double * y2, double * y3, double * y4, double * y5, double * y6, double * y7, double * y8, double * y9, double * y10, int * bidx, int * bsz, double * rhs, int * nrhs, int * mrhs,
     int * nqu, int * nqx, double * rkb, double * hstep, 
     double * Q1, double * Q2, double * R, int * err, double * D, int * iQ2, int * pinfo);

void evdn_(int * dimu, int * dimx, int * stages,
     int * bidx, int * bsz, double * rhs, int * nrhs, int * mrhs,
     int * nqu, int * nqx, double * rkb, double * hstep,
     double * Q1, double * Q2, double * R, int * iQ2, double * DN1, double * DN2, double * DN3, int * lbsz, int * cNidx);

void chlsch_(int * n, int * subn, int * m, int * sb, int * rb, double * A, double * B, int * err,
     int * ns3, int * ns6, int * nbs, int * sdiag, int * ss3, int * ss6, int * sbs);

void evsol_(int * N, int * dimx, int * bidx, int * ynidx, int * bsz, double * rhs, int * nrhs, int * mrhs,
     double * Q1, double * Q2, double * R, int * iQ2, int * P, double * sol);

void evcx_(int * xsize, int * hsize, double * yb, double * sb, double * xi, double * rhs, double * cx1, double * cx2, int * err);

void evpxsz_(int * dimu, double * g,
     int * nlb, int * nub, int * plb, int * pub, int * nr, double * dr,
     double * px, double * ps, double * pz, double * s, double * z, double * mu,
     double * cx1, double * cx2, double * ax1, double * ax2, int * hasB);

void evpxszk_(int * dim, double * px, double * cx1, double * cx2, double * ax1, double * ax2);

void evpy_(int * hsize, double * py, double * cx1, double * cx2, double * ax1, double * ax2, int * hasB);

void merstp_(int * xsize, int * gsize, double * x, double * px, double * s, double * ps, double * alpha, double * xout, double * sout);

void evmer_(int * dim, double * g, double * h,
     int * nlb, int * nub, int * plb, int * pub, int * nr,  double * dr, int * nq,  double * dq,
     double * px, double * ps, double * s, double * mu, double * phi0, double * dphi0);

void evme2_(int * dim, double * g, double * h, int * nlb, int * nub, int * plb, int * pub, int * nr, int * nq, double * s, double * mu, double * phi0);

void evmeru_(int * dimu, int * dimx, int * dimk, int * stages, double * h, double * px, double * dfu, double * dphi0);

void evmerx_(int * dimx, int * dimk, int * stages, double * h, double * px, double * pxA, double * dfx, double * dfxA,
     double * rkb, double * hstep, double * phi0, double * dphi0);

void evmex2_(int * dimx, int * dimk, double * h, double * phi0);

void stepsz_(int * xsize, int * gsize, int * hsize, double * as, double * az, double * dphi0,
     double * px, double * ps, double * pz, double * py, double * tau, double * s, double * z);

void lnsrch_(double * phi0, double * phia, double * dphi0, double * alpha, double * eps, int * res);

void evstu1_(int * dimu, int * dimx, int * dimk,
     double * df, int * nq, double * dq, int * nr, double * dr, int * nlb, int * nub,
     double * x, double * px, double * s, double * ps, double * z, double * pz, double * y1, double * py1, double * y2, double * py2, double * as, double * az,
     double * sk, double * yk, int * err);

void evstx1_(int * dimx, int * dimk,
     double * df, int * nq, double * dq, int * nr, double * dr, int * nlb, int * nub, 
     double * x, double * px, double * s, double * ps, double * z, double * pz, double * y1, double * y2, double * py2, double * as, double * az,
     double * sk, double * yk, int * idf, int * err);

void evstk1_(int * dimx, int * dimk, int * stages,
     double * df, 
     double * x, double * px, double * y, double * as,
     double * sk, double * yk, int * err);

void evstu2_(int * dimu, int * dimx, int * dimk, 
     double * df, int * nq, double * dq, int * nr, double * dr, int * nlb, int * nub,
     double * z, double * y1, double * y2,
     double * sk, double * yk, double * sum);

void evstx2_(int * dimx, int * dimk,
     double * df, int * nq, double * dq, int * nr, double * dr, int * nlb, int * nub,
     double * z, double * y1, double * y2,
     double * sk, double * yk, int * idf, double * sum);

void evstk2_(int * dimx, int * dimk, int * stages,
     double * df,
     double * y,
     double * sk, double * yk, double * sum);

void ifin_(int * gsize, int * hsize, double * g, double * h, double * ieps, double * eps, int * cinf, int * infstp,
     int * stp, int * maxstp, double * mu, int * err, int * iloop);

void ofin_(int * gsize, double * s, double * z, double * mu, double * mtau, double * oeps, double * eps,
     int * stp, int * maxstp, int * ostp, int * oloop);

}
