#ifndef OPTCHOLSCHEDULER_H
#define OPTCHOLSCHEDULER_H

namespace OptChol
{

class Scheduler
{
public:
	virtual ~Scheduler() = 0;
	static int * generate(int n, double * A);

private:

	Scheduler() {}

    typedef struct { int data[2]; } T_CONV_2;
	typedef struct { int data[3]; } T_CONV_3;


};


}


#endif