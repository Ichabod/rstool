#ifndef DISCRGRID_H
#define DISCRGRID_H

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include "vertex.h"
#include "drawable.h"

namespace Discr
{

class Grid
{
public:
	Grid(const int & dim, int * logsize, double * lb, double * ub);
	Grid(const char * filename);
	~Grid();

	inline const int & countNeighbors() const
	{
		return count_neighbors;
	}

	inline const int & countCells() const
	{
		return countall;
	}

	inline unsigned int countMarkedCells() const
	{
		return cells.size();
	}

	inline bool isMarked(const int & cellidx) const
	{
		return mask[cellidx >> 3] & (1 << (cellidx&7));
	}

	inline void setMarked(const int & cellidx)
	{
		if (!isMarked(cellidx))
		{
			mask[cellidx >> 3] |= (1 << (cellidx&7));
			cells.push_back(cellidx);
		}
	}

	inline void setBusy(const int & cellidx, bool isbusy)
	{
		if (isbusy) busy[cellidx >> 3] |= (1 << (cellidx&7));
		else busy[cellidx >> 3] &= (!(unsigned char)(1 << (cellidx&7)));
	}

	inline bool isBusy(const int & cellidx) const
	{
		return busy[cellidx >> 3] & (1 << (cellidx&7));
	}

	inline int getDim() const
	{
		return dim;
	}

	inline double discrStepsize(int dim) const
	{
		return h[dim];
	}

	int getNeighborIdx(const int & cellidx, const int & num);
	int getNeighborIdx(const int & cellidx, const int & num, int * direction);
	void readCellBounds(int cellidx, double * lb, double * ub);
	int getCellIdx(double * x);

	bool save(const char * filename);

	bool generate3DTriangulation(std::vector<Vertex> & vertices, std::vector<Drawable> & triangles,
		bool optimize=true, int smoothdepth=2, double statweight=1.0, double lengthweight=1.2, double cellrestrfactor=1.4151/2.0);



protected:
	void init();
	void generateDrawable(const int & cellidx, int * & corneridx, int * & dimidx, int * & direction, std::vector<Vertex> & vertices, std::map<int, int> & vertices_index, std::vector<Drawable> & drawables);
	void generateDrawables(std::vector<Vertex> & vertices, std::vector<Drawable> & drawables);
	bool isCCW(const Coordinates & p1, const Coordinates & p2, const Coordinates & p3, const Coordinates & normal);
	static void findNeighborVertexes(std::set<int> & found, int vertexnumber, int recursivestep, int maxdepth, std::vector<Vertex> & vertices, std::vector<Drawable> & triangles);


	unsigned char * mask;
	unsigned char * busy;
	double * lb;
	double * ub;
	double * h;
	int * logsize;
	int * size;
	int countall;
	int dim;
	int count_neighbors;
	int * idx_neighbors;
	std::vector<int> cells;

};

}


#endif