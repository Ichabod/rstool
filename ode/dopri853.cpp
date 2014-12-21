#include "dopri853.h"
#include "exception.h"

#include <cstdlib>
#include <cstring>
#include <cfloat>

#define FORTRANNAME(x) x##_
extern "C"
{
	void FORTRANNAME ( dop853 ) ( int * n, ODE::ODEFUNC FCN, double * x, double * y, double * xend, double * rtol, double * atol, int * itol, ODE::SOLOUTFUNC solout, int * iout, double * work, int* lwork, int * iwork, int * liwork, double * rpar, int * ipar, int * did );
}


namespace ODE
{

	DoPri853::DoPri853 ( ODEFUNC func, int ndgl ) : DoPri ( func, ndgl )
	{
		lwork = 11 * ndgl + 21;
		work = ( double* ) malloc ( sizeof ( double ) * lwork );
		if ( work == 0 )
		{
			throw memoryException();
		}

		liwork = 21;
		iwork = ( int* ) malloc ( sizeof ( int ) * liwork );
		if ( work == 0 )
		{
			free ( work );
			throw memoryException();
		}

		DOPRIFUNC = FORTRANNAME ( dop853 );

		setDefaults();

	}

}
