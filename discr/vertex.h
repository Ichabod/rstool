#ifndef DISCRVERTEX_H
#define DISCRVERTEX_H

#include "coordinates.h"

#include <vector>

namespace Discr
{

class Drawable;

class Vertex
{

public:
	Vertex(const int & dim, double * & lb, double * & h, int * & size, const int & idx, const int & number);

	Vertex(const Vertex & in);

	bool operator<( const Vertex & rhs );
	Vertex & operator=( const Vertex & in );

	void readCoordinates(double * x, bool applycorrection=true) const;

	Coordinates coordinates() const;

	inline Coordinates & correction()
	{
		return _correction;
	}

	inline const int & idx() const
	{
		return _idx;
	}

	inline const int & number() const
	{
		return _number;
	}

	inline const int & dim() const
	{
		return _dim;
	}

	void registerDrawable(const int & drawableidx);

	const std::vector<int> & getRegisteredDrawables() const;

private:
	const int _dim;
	double * & _lb;
	double * & _h;
	int * & _size;
	int _idx;
	int _number;
	Coordinates _correction;
	std::vector<int> _usedby;

};

}

#endif