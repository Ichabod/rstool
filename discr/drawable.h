#ifndef DISCRDRAWABLE_H
#define DISCRDRAWABLE_H

#include <cstring>
#include <vector>

#include "vertex.h"
#include "coordinates.h"

namespace Discr
{

class Drawable
{

public:
	Drawable(const std::vector<Vertex> & vertices, const int & count, int * vertexidx, const Coordinates & normal) :
		_vertices(vertices), _count(count), _normal(normal)
	{
		_vertexidx = new int[_count];
		memcpy(_vertexidx, vertexidx, sizeof(int)*_count);
	}

	Drawable(const Drawable & in) : _vertices(in._vertices), _count(in._count), _normal(in._normal)
	{
		_vertexidx = new int[_count];
		memcpy(_vertexidx, in._vertexidx, sizeof(int)*_count);
	}

	Drawable & operator=( const Drawable & in )
	{
		delete[] _vertexidx;
		_count = in._count;

		_vertexidx = new int[_count];
		memcpy(_vertexidx, in._vertexidx, sizeof(int)*_count);

		return *this;
	}

	inline const Coordinates & normal()
	{
		return _normal;
	}

	inline const int & operator[](const int & idx) const
	{
		return _vertexidx[idx];
	}

	inline const int & count() const
	{
		return _count;
	}

private:
	const std::vector<Vertex> & _vertices;
	int _count;
	int * _vertexidx;
	const Coordinates _normal;

};

}

#endif