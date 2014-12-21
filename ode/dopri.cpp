#include "dopri.h"
#include "exception.h"

#include <cfloat>


namespace ODE
{

DoPri::DoPri ( ODEFUNC func, int ndgl ) : Solver ( func, ndgl )
{
}

void DoPri::setDefaults()
{
	work[0] = DBL_EPSILON;
	work[1] = 0.9;
	work[2] = 0.2;
	work[3] = 10.0;
	work[4] = 0.04;
	work[5] = 0.0;
	work[6] = 0.0;

	iwork[0] = 10000000;
	iwork[1] = 1;
	iwork[2] = -1;
	iwork[3] = 1000;
	iwork[4] = 0;
}


void DoPri::calc ( double t, double * rpar, int* ipar )
{

	if ( !initialized ) throw notInitializedException();

	DOPRIFUNC ( &ndgl, f, & ( this->t ), y, &t, &rtol, &atol, &itol, solout, &iout,
	            work, &lwork, iwork, &liwork, rpar, ipar, &did );

	if ( did!=0 )
	{
		switch ( did )
		{
			case -1: throw wrongSolverInputException();
			case -2: throw tooManyStepsException();
			case -3: throw tooSmallStepSizeException();
			case -4: throw stiffException();
		}
	}

	this->t = t;
}

void DoPri::setMaxSteps ( int steps )
{
	if ( steps < 0 ) throw wrongSolverInputException();
	iwork[0] = steps;
}

void DoPri::setStiffnessTest ( int teststeps )
{
	if ( teststeps < 0 ) throw wrongSolverInputException();
	iwork[3] = teststeps;
}

void DoPri::setSafetyFactor ( double safetyfactor )
{
	work[1] = safetyfactor;
}

void DoPri::setStepsizeSelectionParameters ( double fac1, double fac2 )
{
	work[2] = fac1;
	work[3] = fac2;
}

void DoPri::setMaxStepsize ( double stepsize )
{
	if ( stepsize < 0.0 ) throw wrongSolverInputException();
	work[5] = stepsize;
}

void DoPri::setInitialStepsize ( double init )
{
	if ( init < 0.0 ) throw wrongSolverInputException();
	work[6] = init;
}

void DoPri::setBeta ( double beta )
{
	if ( beta < 0.0 || beta > 0.1 ) throw wrongSolverInputException();
	work[4] = beta;
}

int DoPri::getFunctionEvalCount()
{
	return iwork[16];
}

int DoPri::getComputedSteps()
{
	return iwork[17];
}

int DoPri::getAcceptedSteps()
{
	return iwork[18];
}

int DoPri::getRejectedSteps()
{
	return iwork[19];
}

DoPri::~DoPri()
{
}


}


