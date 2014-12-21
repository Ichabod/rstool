#include "initialguess.h"
#include <cstring>

namespace Algorithm
{

int InitialGuess::xsize;
int InitialGuess::gsize;
int InitialGuess::hsize;
int InitialGuess::counter = 0;

InitialGuess::InitialGuess(int count, double xi, double * x, double * z, double * y, double * s_b, double * y_b, double mu)
{
	this->count = count;
	this->x = new double[xsize];
	memcpy(this->x, x, sizeof(double)*xsize);
	if (z != 0)
	{
		this->z = new double[gsize];
		memcpy(this->z, z, sizeof(double)*gsize);
	} else
	{
		this->z = 0;
	}

	if (y != 0)
	{
		this->y = new double[hsize];
		memcpy(this->y, y, sizeof(double)*hsize);
	} else
	{
		this->y = 0;
	}

	if (s_b != 0)
	{
		this->s_b = new double[xsize];
		memcpy(this->s_b, s_b, sizeof(double)*xsize);
	} else
	{
		this->s_b = 0;
	}

	if (y_b != 0)
	{
		this->y_b = new double[xsize];
		memcpy(this->y_b, y_b, sizeof(double)*xsize);
	} else
	{
		this->y_b = 0;
	}

	this->xi = xi;
	this->mu=mu;
	counter++;
}

InitialGuess::InitialGuess()
{
	this->count = count;
	this->x = 0;
	this->z = 0;
	this->y = 0;
	this->s_b = 0;
	this->y_b = 0;

	this->xi = 0;
	this->mu = 0;
	counter++;
}

InitialGuess::~InitialGuess()
{
	if (x != 0) delete[] x;
	if (z != 0) delete[] z;
	if (y != 0) delete[] y;
	if (s_b != 0) delete[] s_b;
	if (y_b != 0) delete[] y_b;
	counter--;
}

void InitialGuess::Destroy(InitialGuess * obj)
{
	if (obj != 0)
	{
		obj->count--;
		if (obj->count <= 0) delete obj;
	}
}

}