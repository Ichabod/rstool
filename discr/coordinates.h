#ifndef DISCRCOORDINATES_H
#define DISCRCOORDINATES_H

#include <cstring>

namespace Discr
{

class Coordinates
{

public:
	Coordinates(const int & dim) : _dim(dim)
	{
		_coordinates = new double[_dim];
	}

	Coordinates(const int & dim, double * pos) : _dim(dim)
	{
		_coordinates = new double[_dim];
		memcpy(_coordinates, pos, sizeof(double)*_dim);
	}

	Coordinates(const int & dim, int * pos) : _dim(dim)
	{
		_coordinates = new double[_dim];
		for (int i=0; i<dim; i++)
		{
			_coordinates[i] = (double)pos[i];
		}
	}

	Coordinates(const Coordinates & in) : _dim(in._dim)
	{
		_coordinates = new double[_dim];
		memcpy(_coordinates, in._coordinates, sizeof(double)*_dim);
	}

	Coordinates & operator=( const Coordinates & in )
	{
		if (_dim == in._dim) memcpy(_coordinates, in._coordinates, sizeof(double)*_dim);
		return *this;
	}

	inline double & operator[](const int & idx) const
	{
		return _coordinates[idx];
	}

	inline const int & dim() const
	{
		return _dim;
	}

private:
	const int _dim;
	double * _coordinates;

};

}

#endif