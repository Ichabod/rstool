#include <iostream>
#include <fstream>
#include <iomanip>
#include <pthread.h>
#include <cstdlib>
#include <chrono>
#include <thread>

#include "discr/grid.h"
#include "algorithm/buffer.h"
#include "deosapp/deos.h"
#include "tools/rtclock.h"

#include "algorithm/reachableset.h"

using namespace std;
using namespace Discr;
using namespace Algorithm;

bool abort_thread=false;

void * statthread(void * arg)
{
	Buffer & buffer = *((Buffer*)arg);
	while (!abort_thread)
	{
		std::this_thread::sleep_for(std::chrono::seconds(2));
		cout << "Buffersize: " << buffer.size() << ", InitialGuess Objects: " << InitialGuess::counter << endl;
	}

	return NULL;
}

void printTimerLine(const char * title, double sec, double totaltime, int markedcells)
{
	cout << setw(13) << title << " : " << fixed << setprecision(6) << sec << " Sek. / " <<
		setprecision(3) << setw(5) << (sec/totaltime)*100.0 << " %, " <<
		setprecision(3) << setw(3) << sec*1000.0/markedcells << " ms/C" << endl;
}

int main()
{
	/****************************
	 * DEOS Problem Parameters  *
	 ****************************/
	int N=30;						// horizon size
	double hstep = 3.0;				// stepsize in seconds

	int logsize[] = { 6,6,6};		// 3D-gridsize (logarithmic, base 2)
	double lb[] = { -2,-5, -3.5};	// grid lower bounds
	double ub[] = { 6,  1.5, 4.2};	// grid upper bounds


	// Target state at time 0
	double targetstart[] = {0, 0, 0, 0, 0, 0,  -0.02, 0.0349, 0.057453, -0.05, 0, 0, 0.99875};

    double c_Mt = 200.0;		// Mass of target
    double c_jt1 = 1000.0;		// inertia torque axis 1
    double c_jt2 = 2000.0;		// inertia torque axis 2
    double c_jt3 = 1000.0;		// inertia torque axis 3
    double c_dtx = 0.0;			// relative x coordinates of docking point
    double c_dty = -1.0;		// relative y coordinates of docking point
    double c_dtz = 0.0;			// relative z coordinates of docking point

	Grid grid(3, logsize, lb,ub);
	Buffer buffer(1);
	DEOSApp::DEOS deos(N,hstep);

    cout << "generating target trajectory...";
    cout.flush();
	if (!deos.generateTarget(targetstart, c_Mt, c_jt1, c_jt2, c_jt3, c_dtx, c_dty, c_dtz))
	{
		cout << "FAILED" << endl;
		return -1;
	}
	cout << "ok" << endl;

	double poshint[] = { 0.0, -1.2, 0.0 };				// guess center of controllability set
	ReachableSet rs(deos, grid, buffer, poshint,8);		// Initialize solver with 8 threads

	pthread_t stat;
	pthread_create(&stat, NULL, statthread, (void *) &buffer);

	Tools::RTClock stopwatch;
	rs.compute();										// Start set computation
	double totalseconds = stopwatch.elapsedSeconds();

	abort_thread=true;
	pthread_join(stat, NULL);

	cout  << endl;
	cout << "TOTAL         : " << fixed << setprecision(6)<< totalseconds << " Sek." << endl;
	cout << "marked cells  : " << grid.countMarkedCells() << endl;
	grid.save("deos_result.grid");
}
