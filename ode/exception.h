#include <exception>

#ifndef ODEEXCEPTION_H
#define ODEEXCEPTION_H

namespace ODE
{
	class Exception : public std::exception
	{
		public:
			virtual void errout ( bool outp = true ) = 0;
	};

	class memoryException : public Exception
	{
		public:
			void errout ( bool outp = true );
	};


	class wrongSolverInputException : public Exception
	{
		public:
			void errout ( bool outp = true );
	};


	class tooManyStepsException : public Exception
	{
		public:
			void errout ( bool outp = true );
	};

	class tooSmallStepSizeException : public Exception
	{
		public:
			void errout ( bool outp = true );
	};


	class notInitializedException : public Exception
	{
		public:
			void errout ( bool outp = true );
	};


	class singularyMatrixException : public Exception
	{
		public:
			void errout ( bool outp = true );
	};

	class stiffException : public Exception
	{
		public:
			void errout ( bool outp = true );
	};

	class unsuccessfulException : public Exception
	{
		public:
			void errout ( bool outp = true );
	};

}

#endif
