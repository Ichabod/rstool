#include "radau5.h"
#include "exception.h"

#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <iostream>

#define FORTRANNAME(x) x##_
extern "C"
{
	void FORTRANNAME ( radau5 ) ( int * n, ODE::ODEFUNC FCN, double * x, double * y, double * xend, double *hinit, double * rtol, double * atol, int * itol, ODE::JACFUNC jac ,int * ijac, int * mljac, int * mujac, ODE::MASSFUNC mas , int * imas, int * mlmas, int * mumas, ODE::SOLOUTFUNC solout, int * iout, double * work, int* lwork, int * iwork, int * liwork, double * rpar, int * ipar, int * did );
}

namespace ODE
{

	Radau5::Radau5 ( ODEFUNC func, int ndgl, JACFUNC jacfunc, int ljac, int ujac, MASSFUNC massfunc, int lmas, int umas ) :
			Radau ( func, ndgl,jacfunc, ljac, ujac, massfunc, lmas, umas )
	{
		int sizejac, sizemas, sizee;

		if ( mljac == ndgl )
		{
			sizejac = ndgl;
			sizee = ndgl;
		}
		else
		{
			sizejac = mljac + mujac + 1;
			sizee = 2*mljac + mujac + 1;
		}

		if ( imas == 0 )
		{
			sizemas = 0;
		}
		else
		{
			if ( mlmas == ndgl )
			{
				sizemas = ndgl;
			}
			else
			{
				sizemas = mlmas + mumas + 1;
			}
		}

		lwork = ndgl* ( sizejac + sizemas + 3*sizee + 12 ) + 20;
		work = ( double* ) malloc ( sizeof ( double ) * lwork );
		if ( work == 0 ) throw memoryException();

		liwork = 3*ndgl + 20;
		iwork = ( int* ) malloc ( sizeof ( int ) * liwork );
		if ( liwork == 0 )
		{
			free ( work );
			throw memoryException();
		}

		RADAUFUNC = FORTRANNAME ( radau5 );

		setDefaults();
	}

	void Radau5::setDefaults()
	{
		iwork[0] = 0;	// hessenberg
		iwork[1] = 10000000; // steps
		iwork[2] = 7; // max newtons
		iwork[3] = 0; // nehmen gute Konvergenz des newtonverfahrens an
		iwork[4] = ndgl;
		iwork[5] = 0;
		iwork[6] = 0;
		iwork[7] = 1; //gustafsson
		iwork[8] = 0;
		iwork[9] = 0;

		work[0] = DBL_EPSILON;
		work[1] = 0.9;
		work[2] = 0.001;
		work[3] = 1E-5;
		work[4] = 1.0;
		work[5] = 1.2;
		work[6] = 0.0;
		work[7] = 0.2;
		work[8] = 8.0;

	}

}
