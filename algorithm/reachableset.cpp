#include <iostream>
#include <iomanip>

#include "reachableset.h"
#include "kernel.h"
#include "../tools/rtclock.h"
#include "nodes/start.h"
#include "nodes/preparesigrhs.h"
#include "nodes/innerblocks.h"
#include "nodes/mergeblocks.h"
#include "nodes/core.h"
#include "nodes/innerloopend.h"
#include "nodes/outerloopend.h"
#include "nodes/finish.h"


using namespace std;
using namespace Algorithm::Nodes;

namespace Algorithm
{

ReachableSet::ReachableSet(Application & app, Discr::Grid & grid, Buffer & buffer, double * const & poshint, const int & threads,
	const double & releps, const int & max_ip_iter, const int & max_divergence_steps) :
	_app(app), _grid(grid), _buffer(buffer), _poshint(poshint),
	_releps(releps), _max_ip_iter(max_ip_iter), _max_divergence_steps(max_divergence_steps)
{

	if (threads<=0)
	{
		_threadcount=1;
	} else
	{
		_threadcount =threads;
	}

	_threadattr = new pthread_attr_t[_threadcount];
	_thread = new pthread_t[_threadcount];
	_threaddone = new bool[_threadcount];
	_storagesize = new int[_threadcount];

	for (int i=0; i<_threadcount; i++)
	{
		_storagesize[i] = 0;
		pthread_attr_init ( &_threadattr[i] );
		pthread_attr_setdetachstate ( &_threadattr[i], PTHREAD_CREATE_JOINABLE );
	}
	pthread_mutex_init(&_condmutex, NULL);
	pthread_cond_init (&_buffercond, NULL);

}

ReachableSet::~ReachableSet()
{
	delete[] _threadattr;
	delete[] _thread;
	delete[] _threaddone;
	pthread_cond_destroy(&_buffercond);
	pthread_mutex_destroy(&_condmutex);

}

void ReachableSet::printTimerLine(const char * title, double sec, double totaltime, int markedcells)
{
	cout << setw(13) << title << " : " << fixed << setprecision(6) << sec << " Sek. / " <<
		setprecision(3) << setw(5) << (sec/totaltime)*100.0 << " %, " <<
		setprecision(3) << setw(3) << sec*1000.0/markedcells << " ms/C" << endl;
}

void * ReachableSet::staticworker(void* ptr)
{
	T_THREADINFO * info = (T_THREADINFO*)(ptr);
	info->base->worker(info->number);

	return 0;
}

void ReachableSet::worker(const int & threadnumber)
{
	_threaddone[threadnumber] = false;

	Kernel kernel(_app.horizonSize(), _app.stepSize(), _app.startTime(), _app.getParameters(), 1e-5,  100, 10);

	Scheduler::Storage storage;

	Start sn(kernel, _buffer, _grid, storage, _threaddone, _buffercond, _condmutex, threadnumber, _threadcount, _storagesize[threadnumber]);
	PrepareSigRhs psrn(kernel);
	InnerBlocks ibn(kernel);
	MergeBlocks mbn(kernel);
	Core cn(kernel);
	InnerLoopEnd iln(kernel);
	OuterLoopEnd oln(kernel);
	Finish fn(kernel, _buffer, _grid, storage, _buffercond, _condmutex, threadnumber);

	sn.setNextNode(&psrn);

	psrn.setNextNode(&ibn);

	ibn.setNextNode(&mbn);
	ibn.setExitNode(&fn);

	mbn.setNextNode(&cn);
	mbn.setExitNode(&fn);

	cn.setNextNode(&iln);
	cn.setExitNode(&fn);

	iln.setExitNode(&fn);
	iln.setLoopNode(&psrn);
	iln.setNextNode(&oln);

	oln.setExitNode(&fn);
	oln.setLoopNode(&psrn);

	fn.setNextNode(&sn);

	Scheduler::Manager manager(sn);

	if (threadnumber == 0)
	{
		InitialGuess * init;

		cout << "generating initial guess...";
		cout.flush();
		if (_app.generateInitialGuess(init, _poshint))
		{
			cout << "FAILED" << endl;
		}
		cout << "ok" << endl;

		init->count = 1;
		init->xi = 2e-1;
		double *startx = init->x+_app.horizonSize()*6;
		cout << "Start = (" << startx[0] << ", " << startx[1] << ", " << startx[2] << ")" << endl;

		_buffer.lock();
		Buffer::ITEM & item = _buffer.push();
		item.guess = init;
		item.grididx = _grid.getCellIdx(startx);
		cout << "Grididx = " << item.grididx << endl;	
		_buffer.unlock();
	}

	Tools::RTClock stopwatch;
	manager.exec();

	double totalseconds = stopwatch.elapsedSeconds();

	_buffer.lock();
	cout << setprecision(5) << endl;
	cout << "THREAD " << threadnumber << ": ELAPSED TIME" << endl;
	cout << "============" << endl << endl;
	printTimerLine("Start", sn.elapsedSeconds(), totalseconds, fn.getCellCounter());
	printTimerLine("PrepareSigRhs", psrn.elapsedSeconds(), totalseconds, fn.getCellCounter());
	printTimerLine("InnerBlocks", ibn.elapsedSeconds(), totalseconds, fn.getCellCounter());
	printTimerLine("MergeBlocks", mbn.elapsedSeconds(), totalseconds, fn.getCellCounter());
	printTimerLine("Core", cn.elapsedSeconds(), totalseconds, fn.getCellCounter());
	printTimerLine("InnerLoopEnd", iln.elapsedSeconds(), totalseconds, fn.getCellCounter());
	printTimerLine("OuterLoopEnd", oln.elapsedSeconds(), totalseconds, fn.getCellCounter());
	printTimerLine("Finish", fn.elapsedSeconds(), totalseconds, fn.getCellCounter());

	_buffer.unlock();

}

void ReachableSet::compute()
{
	T_THREADINFO * info = new T_THREADINFO[_threadcount];

	for (int i=0; i<_threadcount; i++)
	{
		_threaddone[i]=false;
	}
	for (int i=0; i<_threadcount; i++)
	{
		info[i].base=this;
		info[i].number=i;
		pthread_create (&_thread[i], &_threadattr[i], staticworker, &info[i]);
	}

	for (int i=0; i<_threadcount; i++)
	{
		pthread_join(_thread[i],0);
	}

	delete[] info;
}


}