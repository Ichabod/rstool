#include <iostream>
#include "exception.h"

namespace ODE
{
	void memoryException::errout ( bool outp )
	{
		if ( outp )
		{
			std::cerr << "Error occured: memory exception" << std::endl;
		}
	}


	void wrongSolverInputException::errout ( bool outp )
	{
		if ( outp )
		{
			std::cerr << "Error occured: wrong solver input exception" << std::endl;
		}
	}



	void tooManyStepsException::errout ( bool outp )
	{
		if ( outp )
		{
			std::cerr << "Error occured: too many steps" << std::endl;
		}
	}



	void tooSmallStepSizeException::errout ( bool outp )
	{
		if ( outp )
		{
			std::cerr << "Error occured: step size too small" << std::endl;
		}
	}



	void notInitializedException::errout ( bool outp )
	{
		if ( outp )
		{
			std::cerr << "Error occured: parameter not initialized" << std::endl;
		}
	}



	void singularyMatrixException::errout ( bool outp )
	{
		if ( outp )
		{
			std::cerr << "Error occured: matrix singular" << std::endl;
		}
	}


	void stiffException::errout ( bool outp )
	{
		if ( outp )
		{
			std::cerr << "Error occured: function stiff" << std::endl;
		}
	}


	void unsuccessfulException::errout ( bool outp )
	{
		if ( outp )
		{
			std::cerr << "Error occured: routine not successful" << std::endl;
		}
	}
}

