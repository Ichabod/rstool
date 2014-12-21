#ifndef ALGORITHMINITIALGUESS_H
#define ALGORITHMINITIALGUESS_H

namespace Algorithm
{

class InitialGuess
{
public:
	InitialGuess();
	InitialGuess(int count, double xi, double * x, double * z = 0, double * y = 0, double * s_b = 0, double * y_b = 0, double mu = 0.0);
	~InitialGuess();

	static void Destroy(InitialGuess * obj);

	static int xsize;
	static int gsize;
	static int hsize;

	static int counter;

	double xi;
	double * x;
	double * z;
	double * y;
	double * s_b;
	double * y_b;
	double mu;
	int count;

};

}

#endif