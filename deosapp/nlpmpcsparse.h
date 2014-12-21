#ifndef DEOSAPPNLPMPCSparseSPARSE_H
#define DEOSAPPNLPMPCSparseSPARSE_H

#include "IpTNLP.hpp"

using namespace Ipopt;

namespace DEOSApp
{

class NLPMPCSparse : public TNLP
{
public:
	NLPMPCSparse(double hstep, int mpcN, int N, double * target, double * lastx, double * poshint, double * initguess, double c_safety, double * par);

	virtual ~NLPMPCSparse();

	void shift(double * u, double * x, double * k, bool domemshift);



	virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
	                        Index& nnz_h_lag, IndexStyleEnum& index_style);

	virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
	                           Index m, Number* g_l, Number* g_u);

	virtual bool get_starting_point(Index n, bool init_x, Number* x,
	                              bool init_z, Number* z_L, Number* z_U,
	                              Index m, bool init_lambda,
	                              Number* lambda);

	virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

	virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

	virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

	virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
	                      Index m, Index nele_jac, Index* iRow, Index *jCol,
	                      Number* values);

	virtual bool eval_h(Index n, const Number* x, bool new_x,
	                  Number obj_factor, Index m, const Number* lambda,
	                  bool new_lambda, Index nele_hess, Index* iRow,
	                  Index* jCol, Number* values);

	virtual void finalize_solution(SolverReturn status,
	                             Index n, const Number* x, const Number* z_L, const Number* z_U,
	                             Index m, const Number* g, const Number* lambda,
	                             Number obj_value,
	                             const IpoptData* ip_data,
	                             IpoptCalculatedQuantities* ip_cq);

	virtual bool get_scaling_parameters(Number& obj_scaling,
                                    bool& use_x_scaling, Index n,
                                    Number* x_scaling,
                                    bool& use_g_scaling, Index m,
                                    Number* g_scaling);

	private:
		NLPMPCSparse(const NLPMPCSparse&);
		NLPMPCSparse& operator=(const NLPMPCSparse&);

		void recomputeRK(const double * xvec, const double * uvec);
		void rkstepWrapper(double hstep, double * x, double * u, double * par);


		double hstep;
		int mpcN;
		int N;
		double * par;
		double * target;
		double c_safety;
		double * poshint;

		int idx;
		double * lastx;

		double * solution;
		double * z_L;
		double * z_U;
		double * lambda;
		double obj_scaling;
		bool restart;

		int stages;
		double * rkA;
		double * rkb;
		double * rkc;

		double * k;
		double * df;
		int * ip;
		double * dfu;
		double * dfx;
		double * tmpdfu;
		double * tmpdfx;

		double * workfres;
		double * workfres2;
		double * workrktmp;
		double * rkwork;
		int * rkiwork;


};

}


#endif
