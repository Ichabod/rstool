#include "solver.h"

#include <cstdlib>
#include <cstring>

namespace ODE
{

	Solver::Solver ( ODEFUNC func, int ndgl )
	{
		this->ndgl = ndgl;

		f = func;

		itol = 0;
		iout = 0;

		atol = 1E-10;
		rtol = 1E-10;
		solout = 0;

		initialized = false;
		work = 0;
		iwork = 0;

	}

	Solver::~Solver()
	{
		if ( work != 0 ) free ( work );
		if ( iwork != 0 ) free ( iwork );
	}

	void Solver::init ( double t0, double * y0 )
	{
		t = t0;
		y = y0;

		initialized = true;
	}

	void Solver::setTol(double rtol, double atol)
	{
		this->rtol = rtol;
		this->atol = atol;
	}
}

