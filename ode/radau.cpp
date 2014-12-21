#include "radau.h"
#include "exception.h"

#include <cfloat>


namespace ODE
{

Radau::Radau ( ODEFUNC func, int ndgl, JACFUNC jacfunc, int ljac, int ujac, MASSFUNC massfunc, int lmas, int umas ) : Solver ( func, ndgl )
{
	jac = jacfunc;
	mass = massfunc;

	ijac = ( jac == 0 ? 0 : 1 );
	imas = ( mass == 0 ? 0 : 1 );

	mljac = ( ljac == -1 ? ndgl : ljac );
	mujac = ( ujac == -1 ? ndgl : ujac );

	mlmas = ( lmas == -1 ? ndgl : lmas );
	mumas = ( umas == -1 ? ndgl : umas );

	hinit = 0.0;
}


void Radau::calc ( double t, double * rpar, int* ipar )
{

	if ( !initialized ) throw notInitializedException();

	RADAUFUNC ( &ndgl, f, & ( this->t ), y, &t, &hinit, &rtol, &atol, &itol,
	            jac, &ijac, &mljac, &mujac, mass, &imas, &mlmas, &mumas,
	            solout, &iout, work, &lwork, iwork, &liwork, rpar, ipar, &did );

	if ( did!=0 )
	{
		switch ( did )
		{
			case -1: throw wrongSolverInputException();
			case -2: throw tooManyStepsException();
			case -3: throw tooSmallStepSizeException();
			case -4: throw singularyMatrixException();
		}
	}

	this->t = t;
}

void Radau::setMaxSteps ( int steps )
{
	if ( steps < 0 ) throw wrongSolverInputException();
	iwork[1] = steps;
}

void Radau::setSafetyFactor ( double safetyfactor )
{
	work[1] = safetyfactor;
}

void Radau::setStepsizeSelectionParameters ( double fac1, double fac2 )
{
	work[7] = fac1;
	work[8] = fac2;
}

void Radau::setMaxStepsize ( double stepsize )
{
	if ( stepsize < 0.0 ) throw wrongSolverInputException();
	work[6] = stepsize;
}

void Radau::setInitialStepsize ( double init )
{
	if ( init < 0.0 ) throw wrongSolverInputException();
	hinit = init;
}

int Radau::getFunctionEvalCount()
{
	return iwork[13];
}

int Radau::getComputedSteps()
{
	return iwork[15];
}

int Radau::getAcceptedSteps()
{
	return iwork[16];
}

int Radau::getRejectedSteps()
{
	return iwork[17];
}

Radau::~ Radau()
{
}

void Radau::enableHessenbergTransformation ( bool enable )
{
	iwork[0] = ( enable == true? 1 : 0 );
}

void Radau::setMaxNewtonSteps ( int steps )
{
	if ( steps < 0 ) throw wrongSolverInputException();
	iwork[2] = steps;
}

void Radau::enableStartingValues ( bool enable )
{
	iwork[3] = ( enable == true? 1 : 0 );
}

bool Radau::badConvergence()
{
	return getComputedSteps() > getAcceptedSteps() + getRejectedSteps();
}

void Radau::setAlgebraicSystemParameters ( int dim_index1, int dim_index2, int dim_index3, int m1, int m2 )
{
	iwork[4] = dim_index1;
	iwork[5] = dim_index2;
	iwork[6] = dim_index3;
	iwork[8] = m1;
	iwork[9] = m2;
}

void Radau::setJacobiRecomputation ( double step )
{
	work[2] = step;
}

void Radau::setStepsizeStrategy ( STEPSIZESTRATEGY strategy )
{
	iwork[7] = ( strategy == Gustaffson ? 1 : 2 );
}

void Radau::setStepsizeFixingParameters ( double fac1, double fac2 )
{
	work[4] = fac1;
	work[5] = fac2;
}


}
