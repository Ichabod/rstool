#include "grid.h"
#include "vertex.h"
#include "gridoptnlp.h"

#include <IpIpoptApplication.hpp>

#include <cmath>
#include <iostream>
#include <cstring>
#include <fstream>
#include <map>
#include <set>

using namespace std;

namespace Discr
{

Grid::Grid(const int & dim, int * logsize, double * lb, double * ub)
{
	this->dim = dim;
	if (dim < 1 || dim > 4) throw 1;

	this->logsize = new int[dim];
	this->lb = new double[dim];
	this->ub = new double[dim];

	for (int i=0; i<dim; i++)
	{
		this->logsize[i] = logsize[i];
		this->lb[i] = lb[i];
		this->ub[i] = ub[i];
	}

	init();
}

Grid::Grid(const char * filename)
{
	ifstream file(filename);
	if (!file.good() || file.eof()) throw 1;

	unsigned char fileversion;
	file.read((char*)&fileversion, 1);
	if (fileversion != 1) throw 1;

	file.read((char*)&dim, sizeof(int));
	if (dim < 1 || dim > 4) throw 1;

	this->logsize = new int[dim];
	this->lb = new double[dim];
	this->ub = new double[dim];

	file.read((char*)logsize, sizeof(int)*dim);
	file.read((char*)lb, sizeof(double)*dim);
	file.read((char*)ub, sizeof(double)*dim);

	int size;
	file.read((char*)&size, sizeof(int));
	cells.resize(size);

	file.read((char*)&cells[0], sizeof(int)*size);
	if (file.fail() || file.bad()) throw 1;

	file.close();

	init();

	for (int i=0; i<size; i++)
	{
		mask[cells[i] >> 3] |= (1 << (cells[i]&7));
	}
}

Grid::~Grid()
{
	delete[] logsize;
	delete[] size;
	delete[] lb;
	delete[] ub;
	delete[] h;
	delete[] idx_neighbors;
	delete[] mask;
}

void Grid::init()
{
	this->size = new int[dim];
	this->h = new double[dim];
	countall = 1;

	for (int i=0; i<dim; i++)
	{
		size[i] = (1<<logsize[i]);
		h[i] = (ub[i] - lb[i]) / (double)size[i];
		countall *= size[i];
	}

	count_neighbors = round(pow(3.0, dim))-1;
	idx_neighbors = new int[count_neighbors];
	for (int i=1; i<=count_neighbors; i++)
	{
		int relative=0;
		int tmp = i;
		int ofs=1;
		for (int j=0; j<dim; j++)
		{
			int d = tmp%3;
			if (d == 2) d = -1;
			relative += d*ofs;
			tmp /=3;

			ofs *= size[j];

		}
		idx_neighbors[i-1] = relative;
	}
	mask = new unsigned char[(countall>>3)+1];
	busy = new unsigned char[(countall>>3)+1];
	memset(mask, 0, sizeof(unsigned char)*((countall>>3)+1));
	memset(busy, 0, sizeof(unsigned char)*((countall>>3)+1));
}


void Grid::readCellBounds(int cellidx, double * lb, double * ub)
{
	for (int i=0; i<dim; i++)
	{
		int dimidx = cellidx & ((1<<logsize[i])-1);
		cellidx >>= logsize[i];
		lb[i] = this->lb[i] + h[i]*(double)dimidx;
		ub[i] = this->lb[i] + h[i]*(double)(dimidx+1);
	}
}

int Grid::getCellIdx(double * x)
{
	int res = 0;
	int ofs = 1;
	for (int i=0; i<dim; i++)
	{
		int dimidx = floor((x[i] - lb[i])/h[i]);
		if (dimidx < 0 || dimidx >= (int)size[i])
		{
			return -1;
		}
		res += dimidx*ofs;
		ofs *= size[i];
	}
	return res;
}

int Grid::getNeighborIdx(const int & cellidx, const int & num, int * direction)
{
	int res = cellidx+idx_neighbors[num];
	if (res < 0 || res >= countall) return -1;
	int logofs = 0;
	for (int i=0; i<dim; i++)
	{
		int bitmask = ((1<<logsize[i])-1)<<logofs;
		int diff = ((res & bitmask)>>logofs) - ((cellidx & bitmask)>>logofs);
		logofs += logsize[i];
		if ( diff < -1 || diff > 1)
		{
			return -1;
		}
		if (direction)
		{
			direction[i] = diff;
		}
	}

	return res;
}

int Grid::getNeighborIdx(const int & cellidx, const int & num)
{
	return getNeighborIdx(cellidx,num,0);
}

bool Grid::save(const char * filename)
{
	ofstream file(filename);
	if (!file.good()) return false;

	unsigned char fileversion = 1;

	file.write((const char*)&fileversion, 1);
	file.write((const char*)&dim, sizeof(int));
	file.write((const char*)logsize, sizeof(int)*dim);
	file.write((const char*)lb, sizeof(double)*dim);
	file.write((const char*)ub, sizeof(double)*dim);
	int size = cells.size();
	file.write((const char*)&size, sizeof(int));
	file.write((const char*)&cells[0], sizeof(int)*size);

	file.close();

	return true;
}

void Grid::generateDrawable(const int & cellidx, int * & corneridx, int * & dimidx, int * & direction, vector<Vertex> & vertices, map<int, int> & vertices_index, vector<Drawable> & drawables )
{
//	cout << "dir=" << direction[0] << " " << direction[1] << " " << direction[2] << " " << endl;

	// Knoten ist nicht diagonal benachbart?
	int dircounter=0;
	for (int i=0; i<dim; i++)
	{
		if (direction[i] != 0) dircounter++;
	}
	if (dircounter == 1)
	{

		// Alle Ecken durchlaufen
		int idxcounter=0;
		for (int i=0; i<(1<<dim); i++)
		{
			bool validnode=true;
			// Testen, ob die Ecke in beiden Zellen vorkommt
			for (int j=0; j<dim && validnode; j++)
			{
				if (  abs( ((i&(1<<j))!=0?1:-1) - direction[j])>1  )
				{
					validnode=false;
				}
			}

			if (validnode)
			{
				int idx = 0;
				int shift=1;
				for (int j=0; j<dim; j++)
				{
					idx += shift*((i&(1<<j))==0 ? dimidx[j] : dimidx[j]+1);
					shift *= size[j]+1;
				}

				map<int,int>::iterator it = vertices_index.find(idx);
				if (it == vertices_index.end())
				{
					int arrpos = vertices.size();
					pair<map<int,int>::iterator,bool> ret = vertices_index.insert(pair<int,int>(idx, arrpos));
					it = ret.first;
					vertices.push_back(Vertex(dim, lb, h, size, idx, arrpos));
				}

				corneridx[idxcounter++] = it->second;
			}
		}

		drawables.push_back(Drawable(vertices, idxcounter, corneridx, Coordinates(dim, direction)));
	}

}

void Grid::generateDrawables(vector<Vertex> & vertices, vector<Drawable> & drawables)
{
	int * corneridx = new int[1<<(dim-1)];
	int * dimidx = new int[dim];
	int * direction = new int[dim];
	map<int, int> vertices_index;

	for (unsigned int ci=0; ci<cells.size(); ci++)
	{
		int cellidx = cells[ci];
		int tempidx = cellidx;
		for (int j=0; j<dim; j++)
		{
			dimidx[j] = tempidx & ((1<<logsize[j])-1);
			tempidx >>= logsize[j];
		}

		// Für alle Nachbarn
		for (int nbi=0; nbi<countNeighbors(); nbi++)
		{
			int nbidx = getNeighborIdx(cellidx, nbi, direction);
			if (nbidx != -1 && !isMarked(nbidx))
			{
				// Nachbar ist leer, Knoten sind möglich
				generateDrawable(cellidx, corneridx, dimidx, direction, vertices, vertices_index, drawables );
			}
		}

		// Und am Rand
		for (int i=0; i<dim; i++)
		{
			for (int j=0; j<dim; j++) direction[j] = 0;

			bool isborder = false;
			if (dimidx[i] == 0)
			{
				direction[i] = -1;
				isborder = true;
			}
			else if (dimidx[i] == size[i]-1)
			{
				direction[i] = 1;
				isborder = true;
			}

			if (isborder) generateDrawable(cellidx, corneridx, dimidx, direction, vertices, vertices_index, drawables );
		}
	}
	delete[] direction;
	delete[] dimidx;
	delete[] corneridx;

}

bool Grid::generate3DTriangulation(vector<Vertex> & vertices, vector<Drawable> & triangles, bool optimize, int smoothdepth, double statweight, double lengthweight, double cellrestrfactor)
{
	if (dim != 3) return false;

	vertices.clear();
	triangles.clear();

	vector<Drawable> drawables;
	
	Coordinates points[] = { Coordinates(3), Coordinates(3), Coordinates(3), Coordinates(3) };

	generateDrawables(vertices, drawables);

	for (unsigned int i=0; i<drawables.size(); i++)
	{
		if (drawables[i].count() != 4)
		{
			cout << "drawables must have 4 vertices... something went wrong..." << endl;
			return false;
		}

		const Coordinates & normal = drawables[i].normal();

		for (int j=0; j<4; j++)	{
			points[j] = vertices[drawables[i][j]].coordinates();
/*			cout << "( ";
			for (int j1=0; j1<3; j1++)
			{
				cout << points[j][j1] << " ";
			}
			cout << ")" << endl;*/
		}
//		cout << endl;

		bool found=false;
		int j1,j2,j3,j4;
		for (j1=0; j1<4; j1++)
		{
			for (j2=0; j2<4; j2++)
			{
				if (j1 != j2)
				{
					for (j3=0; j3<4; j3++)
					{
						if (j1 != j3 && j2 != j3 && isCCW(points[j1], points[j2], points[j3], normal))
						{
							for (j4=0; j4<4; j4++)
							{
								if (j1 != j4 && j2 != j4 && j3 != j4)
								{
									if (isCCW(points[j1], points[j4], points[j2], normal))
									{
										found = true;
										break;
									}
								}
							}
							if (found) break;
						}
					}
					if (found) break;
				}
			}
			if (found) break;
		}

		if (found)
		{
//			cout << j1 << " " << j2 << " " << j3 << " " << j4 << endl;
			int idx1[] = { drawables[i][j1], drawables[i][j2], drawables[i][j3]};
			triangles.push_back(Drawable(vertices, 3, idx1, normal));

			vertices[idx1[0]].registerDrawable(triangles.size()-1);
			vertices[idx1[1]].registerDrawable(triangles.size()-1);
			vertices[idx1[2]].registerDrawable(triangles.size()-1);

			int idx2[] = { drawables[i][j1], drawables[i][j4], drawables[i][j2]};
			triangles.push_back(Drawable(vertices, 3, idx2, normal));

			vertices[idx2[0]].registerDrawable(triangles.size()-1);
			vertices[idx2[1]].registerDrawable(triangles.size()-1);
			vertices[idx2[2]].registerDrawable(triangles.size()-1);
		}
	}

	if (smoothdepth > 0)
	{
		// smoothen the mesh...
		double pv1[3];
		double pv2[3];
		set<int> found;
		for (unsigned int i=0; i<vertices.size(); i++)
		{
			found.clear();
			found.insert(i);
			findNeighborVertexes(found, i, 0, smoothdepth, vertices, triangles);
			Vertex & v1 = vertices[i];
			double * corr = &(v1.correction()[0]);
			for (int l=0; l<3; l++)
			{
				corr[l] = 0.0;
			}
			v1.readCoordinates(pv1, false);

		    for (set<int>::iterator it = found.begin(); it != found.end(); it++)
		    {
				Vertex & v2 = vertices[*it];
				if (v1.number() != v2.number())
				{
					v2.readCoordinates(pv2, false);
					for (int l=0; l<3; l++)
					{
						corr[l] += pv2[l] - pv1[l];
					}
				}
		    }
			for (int l=0; l<3; l++)
			{
				corr[l] /= (double)(found.size());
			}
		}
	}

	if (optimize)
	{
		Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptapp = IpoptApplicationFactory();
	    ipoptapp->Initialize();

	    ipoptapp->Options()->SetStringValue ( "sb", "yes" );
	    ipoptapp->Options()->SetIntegerValue ( "print_level", 5 );
	    ipoptapp->Options()->SetStringValue ( "linear_solver", "ma57" );
	//    ipoptapp->Options()->SetStringValue ( "hessian_approximation", "limited-memory" );
//	    ipoptapp->Options()->SetStringValue ( "derivative_test", "second-order");
	    ipoptapp->Options()->SetStringValue ( "hessian_constant", "yes");
	    ipoptapp->Options()->SetStringValue ("nlp_scaling_method", "user-scaling");
	    ipoptapp->Options()->SetNumericValue ( "tol", 1E-7 );
	   
	    SmartPtr<GridOptNLP> nlp = new GridOptNLP(*this, vertices, triangles, statweight, lengthweight, cellrestrfactor);

	    ApplicationReturnStatus  status = ipoptapp->OptimizeTNLP(nlp);
	    if (status < 0)
	    {
	    	cout << "FAILED" << endl;
	        return false;
	    }


	}


	return true;
}

void Grid::findNeighborVertexes(set<int> & found, int vertexnumber, int recursivestep, int maxdepth, vector<Vertex> & vertices, vector<Drawable> & triangles)
{
	recursivestep++;
	Vertex & v1 = vertices[vertexnumber];
	for (unsigned int j=0; j<v1.getRegisteredDrawables().size(); j++)
	{
	    int drawableidx = v1.getRegisteredDrawables()[j];

	    for (int k=0; k<3; k++)
	    {
	    	Vertex & v2 = vertices[triangles[drawableidx][k]];
	    	if (v1.number() != v2.number())
	    	{
	    		found.insert(v2.number());
	    		if (recursivestep < maxdepth)
	    		{
	    			findNeighborVertexes(found, v2.number(), recursivestep, maxdepth, vertices, triangles);
	    		}
	    	}
	    }
	}	
}

bool Grid::isCCW(const Coordinates & p1, const Coordinates & p2, const Coordinates & p3, const Coordinates & normal)
{
	double cross[3];
	double v1[] = 
	{
		p2[0]-p1[0],
		p2[1]-p1[1],
		p2[2]-p1[2]
	};
	double v2[] =
	{
		p3[0]-p2[0],
		p3[1]-p2[1],
		p3[2]-p2[2]
	};
	for (int j=0; j<3; j++)
	{
		cross[j] = v1[(j+1)%3]*v2[(j+2)%3] - v1[(j+2)%3]*v2[(j+1)%3];
	}
	return cross[0]*normal[0] + cross[1]*normal[1] + cross[2]*normal[2] > 0.0;

}


}