#ifndef DISCRGRIDOPTNLP_H
#define DISCRGRIDOPTNLP_H

#include <vector>
#include <map>

#include "IpTNLP.hpp"

#include "grid.h"
#include "vertex.h"
#include "drawable.h"

using namespace Ipopt;

namespace Discr
{

class GridOptNLP  : public TNLP
{
public:
	GridOptNLP (const Grid & grid, std::vector<Vertex> & vertices, const std::vector<Drawable> & triangles, const double & statweight, const double & lengthweight, const double & cellrestrfactor);

	virtual ~GridOptNLP ();

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
		GridOptNLP (const GridOptNLP &);
		GridOptNLP & operator=(const GridOptNLP &);
		void hesshelper_layout(int row, int col, int & cursor, double val);
		void hesshelper_generate();

		const Grid & _grid;
		std::vector<Vertex> & _vertices;
		const std::vector<Drawable> & _triangles;
		double _discrweight[3];
		const double _statweight;
		const double _lengthweight;
		const double _cellrestrfactor;
		double _mindiscr;

		typedef struct
		{
			int cursor;
			double value;
		} T_HESS;
		std::map<int,T_HESS> _hessmap;
};

}


#endif
