#include "vertex.h"
#include "drawable.h"

namespace Discr
{

Vertex::Vertex(const int & dim, double * & lb, double * & h, int * & size, const int & idx, const int & number) :
	_dim(dim), _lb(lb), _h(h), _size(size), _idx(idx), _number(number), _correction(dim)
{
	for (int i=0; i<dim; i++)
	{
		_correction[i] = 0.0;
	}
}

Vertex::Vertex(const Vertex & in) : _dim(in._dim), _lb(in._lb), _h(in._h), _size(in._size), _idx(in._idx) , _number(in._number), _correction(in._correction) {}

bool Vertex::operator<( const Vertex & rhs ) { return _idx < rhs._idx; }

Vertex & Vertex::operator=( const Vertex & in )
{
	_idx = in._idx;
	_number = in._number;
	return *this;
}

void Vertex::readCoordinates(double * x, bool applycorrection) const
{
	int tmpidx = _idx;
	for (int i=0; i<_dim; i++)
	{
		int dimidx = tmpidx  % (_size[i]+1);
		tmpidx /= _size[i]+1;
		x[i] = _lb[i] + _h[i]*(double)dimidx;
		if (applycorrection) x[i] += _correction[i];
	}
}

Coordinates Vertex::coordinates() const
{
	Coordinates res(_dim);
	readCoordinates(&res[0]);

	return res;
}

void Vertex::registerDrawable(const int & drawableidx)
{
	_usedby.push_back(drawableidx);
}

const std::vector<int> & Vertex::getRegisteredDrawables() const
{
	return _usedby;
}

}

