#ifndef ALGORITHM_APPLICATION_H
#define ALGORITHM_APPLICATION_H

#include "initialguess.h"

namespace Algorithm
{

class Application
{
public:

	virtual ~Application() {}
	virtual int generateInitialGuess(Algorithm::InitialGuess * & initialguess, double * poshint) = 0;
	virtual int horizonSize() const = 0;
	virtual int stepSize() const = 0;
	virtual double startTime() const = 0;
	virtual double * getParameters() = 0;

//	virtual int dimU() const = 0;
//	virtual int dimX() const = 0;
//	virtual int stages() const = 0;
//	virtual void rungeKuttaMethod(double * rkA, double * rkb, double * rkc) const = 0;

};


}

#endif
