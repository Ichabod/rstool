#ifndef ALGORITHM_REACHABLESET_H
#define ALGORITHM_REACHABLESET_H

#include "application.h"
#include "buffer.h"
#include "../discr/grid.h"


namespace Algorithm
{

class ReachableSet
{

public:
	ReachableSet(Application & app, Discr::Grid & grid, Buffer & buffer, double * const & poshint, const int & threads=0,
		const double & releps = 1e-2, const int & max_ip_iter = 100, const int & max_divergence_steps = 10);

	~ReachableSet();

	void compute();

protected:
	typedef struct
	{
		ReachableSet * base;
		int number;
	} T_THREADINFO;

	void printTimerLine(const char * title, double sec, double totaltime, int markedcells);

	static void * staticworker(void* ptr);
	void worker(const int & threadnumber);

	Application & _app;
	Discr::Grid & _grid;
	Buffer & _buffer;
	int _threadcount;
	double * const & _poshint;
	const double & _releps;
	const int & _max_ip_iter;
	const int & _max_divergence_steps;

	pthread_attr_t * _threadattr;
	pthread_t * _thread;
	bool * _threaddone;
	int * _storagesize;
	pthread_cond_t _buffercond;
	pthread_mutex_t _condmutex;



};


}


#endif