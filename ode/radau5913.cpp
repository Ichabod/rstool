#include "radau5913.h"
#include "exception.h"

#include <cstdlib>
#include <cstring>
#include <cfloat>

#define FORTRANNAME(x) x##_
extern "C"
{
	void FORTRANNAME ( radau ) ( int * n, ODE::ODEFUNC FCN, double * x, double * y, double * xend, double *hinit, double * rtol, double * atol, int * itol, ODE::JACFUNC jac ,int * ijac, int * mljac, int * mujac, ODE::MASSFUNC mas , int * imas, int * mlmas, int * mumas, ODE::SOLOUTFUNC solout, int * iout, double * work, int* lwork, int * iwork, int * liwork, double * rpar, int * ipar, int * did );
}

namespace ODE
{

Radau5913::Radau5913 ( ODEFUNC func, int ndgl, JACFUNC jacfunc, int ljac, int ujac, MASSFUNC massfunc, int lmas, int umas ) :
		Radau ( func, ndgl,jacfunc, ljac, ujac, massfunc, lmas, umas )
{
	maxstages = 7;

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

	lwork = ndgl* ( sizejac + sizemas + maxstages*sizee + 3*maxstages + 3 ) + 20;
	work = ( double* ) malloc ( sizeof ( double ) * lwork );
	if ( work == 0 ) throw memoryException();

	liwork = ( 2 + ( maxstages-1 ) /2 ) * ndgl + 20;
	iwork = ( int* ) malloc ( sizeof ( int ) * liwork );
	if ( liwork == 0 )
	{
		free ( work );
		throw memoryException();
	}

	RADAUFUNC = FORTRANNAME ( radau );

	setDefaults();
}

void Radau5913::setDefaults()
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
	iwork[10] = 3;
	iwork[11] = maxstages;
	iwork[12] = 3;

	work[0] = DBL_EPSILON;
	work[1] = 0.9;
	work[2] = 0.001;
	work[3] = 1E-5;
	work[4] = 1.0;
	work[5] = 1.2;
	work[6] = 0.0;
	work[7] = 0.2;
	work[8] = 8.0;
	work[9] = 0.002;
	work[10] = 0.8;
	work[11] = 1.2;
	work[12] = 0.8;

}

void Radau5913::setOrderSelectionParameters ( double inc_factor, double dec_factor, double hfac1, double hfac2 )
{
	work[9] = inc_factor;
	work[10] = dec_factor;
	work[11] = hfac1;
	work[12] = hfac2;
}

void Radau5913::setStages ( int minstage, int maxstage, int firststage )
{
	if ( !isValidStage ( minstage ) || !isValidStage ( maxstage ) || !isValidStage ( firststage ) )
	{
		throw wrongSolverInputException();
	}

	iwork[10] = minstage;
	iwork[11] = maxstage;
	iwork[12] = firststage;
}

bool Radau5913::isValidStage ( int stage )
{
	return ( stage == 1 || stage == 3 || stage == 5 || stage == 7 );
}

}

