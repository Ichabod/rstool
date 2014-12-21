#include "gridoptnlp.h"
#include <cstring>
#include <cmath>
#include <cassert>
#include <map>

using namespace Ipopt;
using namespace std;

namespace Discr
{

GridOptNLP::GridOptNLP(const Grid & grid, vector<Vertex> & vertices, const vector<Drawable> & triangles, const double & statweight, const double & lengthweight, const double & cellrestrfactor) :
    _grid(grid), _vertices(vertices), _triangles(triangles), _statweight(statweight), _lengthweight(lengthweight), _cellrestrfactor(cellrestrfactor)
{
    _mindiscr = min(grid.discrStepsize(0), min(grid.discrStepsize(1), grid.discrStepsize(2)));

    for (int i=0; i<3; i++)
    {
//        _discrweight[i] = pow(_mindiscr/grid.discrStepsize(i), 2.0);
        _discrweight[i] = pow(1.0/grid.discrStepsize(i), 2.0);
    }
}

GridOptNLP::~GridOptNLP()
{
}

bool GridOptNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    n = _vertices.size()*3;
    m = 0;

    nnz_jac_g = 0;

    hesshelper_generate();

    nnz_h_lag = 3*_vertices.size() + _hessmap.size();

    index_style = FORTRAN_STYLE;

    return true;
}

bool GridOptNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
    for (unsigned int i=0; i<_vertices.size(); i++)
    {
        for (int j=0; j<3; j++)
        {
            x_l[i*3+j] = -_grid.discrStepsize(j)*_cellrestrfactor;
            x_u[i*3+j] = _grid.discrStepsize(j)*_cellrestrfactor;
        }
    }

    return true;
}

bool GridOptNLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
    for (unsigned int i=0; i<_vertices.size()*3; i++)
    {
        x[i] = 0.0;
    }

    return true;
}

bool GridOptNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    double d1[] = {0.0, 0.0, 0.0};

    for (unsigned int i=0; i<_vertices.size(); i++)
    {
        d1[0] += x[i*3  ]*x[i*3  ];
        d1[1] += x[i*3+1]*x[i*3+1];
        d1[2] += x[i*3+2]*x[i*3+2];
    }

    double d2 = 0.0;
    for (unsigned int i=0; i<_triangles.size(); i++)
    {
        double p1[3];
        double p2[3];
        double p3[3];

        _vertices[_triangles[i][0]].readCoordinates(p1);
        _vertices[_triangles[i][1]].readCoordinates(p2);
        _vertices[_triangles[i][2]].readCoordinates(p3);

        for (int j=0; j<3; j++)
        {
            p1[j] += x[_triangles[i][0]*3+j];
            p2[j] += x[_triangles[i][1]*3+j];
            p3[j] += x[_triangles[i][2]*3+j];
        }

        for (int j=0; j<3; j++)
        {
            d2 += (p1[j]-p2[j])*(p1[j]-p2[j])*_discrweight[j];
            d2 += (p1[j]-p3[j])*(p1[j]-p3[j])*_discrweight[j];
            d2 += (p2[j]-p3[j])*(p2[j]-p3[j])*_discrweight[j];
        }
    }

    obj_value = _statweight* ( _discrweight[0]*d1[0] + _discrweight[1]*d1[1] + _discrweight[2]*d1[2] ) + _lengthweight*d2;
    return true;
}

bool GridOptNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    for (unsigned int i=0; i<_vertices.size(); i++)
    {
        for (int j=0; j<3; j++)
        {
            grad_f[i*3+j] = _statweight*_discrweight[j]*2.0*x[i*3+j];
        }
    }

    for (unsigned int i=0; i<_triangles.size(); i++)
    {
        double p1[3];
        double p2[3];
        double p3[3];

        _vertices[_triangles[i][0]].readCoordinates(p1);
        _vertices[_triangles[i][1]].readCoordinates(p2);
        _vertices[_triangles[i][2]].readCoordinates(p3);

        for (int j=0; j<3; j++)
        {
            p1[j] += x[_triangles[i][0]*3+j];
            p2[j] += x[_triangles[i][1]*3+j];
            p3[j] += x[_triangles[i][2]*3+j];
        }

        for (int j=0; j<3; j++)
        {
            grad_f[_triangles[i][0]*3+j] += _lengthweight* (  2.0*(p1[j]-p2[j]) + 2.0*(p1[j]-p3[j]) )*_discrweight[j];
            grad_f[_triangles[i][1]*3+j] += _lengthweight* ( -2.0*(p1[j]-p2[j]) + 2.0*(p2[j]-p3[j]) )*_discrweight[j];
            grad_f[_triangles[i][2]*3+j] += _lengthweight* ( -2.0*(p1[j]-p3[j]) - 2.0*(p2[j]-p3[j]) )*_discrweight[j];
        }
    }
    return true;
}

bool GridOptNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    return false;
}

bool GridOptNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
    return false;
}

void GridOptNLP::hesshelper_layout(int row, int col, int & cursor, double val)
{
    int _row,_col;

    // lower tri
    _row = max(row , col);
    _col = min(row , col);

    int key = (_row<<16) + _col;
    map<int, T_HESS>::iterator it = _hessmap.find(key);
    if (it == _hessmap.end())
    {
        T_HESS item;
        item.cursor = cursor;
        item.value = _lengthweight*val;
        _hessmap.insert(pair<int,T_HESS>(key, item));
        cursor++;
    } else
    {
        it->second.value += _lengthweight* val;
    }
}

void GridOptNLP::hesshelper_generate()
{
    _hessmap.clear();
    int cursor = _vertices.size()*3;
    for (unsigned int i=0; i<_triangles.size(); i++)
    {
        for (int j=0; j<3; j++)
        {
            hesshelper_layout(_triangles[i][0]*3+j, _triangles[i][1]*3+j, cursor, -2.0*_discrweight[j]);
            hesshelper_layout(_triangles[i][0]*3+j, _triangles[i][2]*3+j, cursor, -2.0*_discrweight[j]);
            hesshelper_layout(_triangles[i][1]*3+j, _triangles[i][2]*3+j, cursor, -2.0*_discrweight[j]);
        }
    }
}

bool GridOptNLP::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
    if (values == 0)
    {
        int cursor = 0;
        for (unsigned int i=0; i<_vertices.size(); i++)
        {
            for (int j=0; j<3; j++)
            {
                if (iRow != 0 && jCol != 0)
                {
                    iRow[cursor] = jCol[cursor]= i*3+j+1; // includes diag from triangles
                }
                cursor++;
            }
        }

        for (map<int,T_HESS>::iterator it = _hessmap.begin(); it != _hessmap.end(); it++)
        {
            int col = it->first & 0x0000FFFF;
            int row = (it->first & 0xFFFF0000) >> 16;
            iRow[it->second.cursor] = row+1;
            jCol[it->second.cursor] = col+1;

            cursor++;
        }
        assert(cursor == nele_hess);
    } else
    {
        int cursor = 0;
        for (unsigned int i=0; i<_vertices.size(); i++)
        {
            for (int j=0; j<3; j++)
            {
                values[cursor++] = obj_factor*_statweight*_discrweight[j]*2.0;
            }
        }

        for (int i=cursor; i<nele_hess; i++)
        {
            values[i] = 0.0;
        }

        for (unsigned int i=0; i<_triangles.size(); i++)
        {
            for (int j=0; j<3; j++)
            {
                values[_triangles[i][0]*3+j] += obj_factor*_lengthweight* ( 2.0 + 2.0 )*_discrweight[j];
                values[_triangles[i][1]*3+j] += obj_factor*_lengthweight* ( 2.0 + 2.0 )*_discrweight[j];
                values[_triangles[i][2]*3+j] += obj_factor*_lengthweight* ( 2.0 + 2.0 )*_discrweight[j];
            }
        }
        for (map<int,T_HESS>::iterator it = _hessmap.begin(); it != _hessmap.end(); it++)
        {
            values[it->second.cursor] = obj_factor*it->second.value;
        }
    }


    return true;
}

void GridOptNLP::finalize_solution(SolverReturn status,
                            Index n, const Number* x, const Number* z_L, const Number* z_U,
                            Index m, const Number* g, const Number* lambda,
                            Number obj_value,
                            const IpoptData* ip_data,
                            IpoptCalculatedQuantities* ip_cq)
{
    for (unsigned int i=0; i<_vertices.size(); i++)
    {
        Coordinates & correction = _vertices[i].correction();
        for (int j=0; j<3;j++)
        {
            correction[j] += x[i*3+j];
        }
    }
}

bool GridOptNLP::get_scaling_parameters(Number& obj_scaling,
                                    bool& use_x_scaling, Index n,
                                    Number* x_scaling,
                                    bool& use_g_scaling, Index m,
                                    Number* g_scaling)
{
    use_x_scaling=false;
    use_g_scaling=false;

    obj_scaling=1.0;
//    obj_scaling=1.0/(_mindiscr*_mindiscr) * 1e-3;

    return true;
}

}

