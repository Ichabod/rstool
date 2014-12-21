#include "outerloopend.h"

#include <iostream>
#include <cmath>

using namespace std;
using namespace Scheduler;

//#define VERBOSE

namespace Algorithm { namespace Nodes
{


OuterLoopEnd::OuterLoopEnd(Algorithm::Kernel & kernel) :
	Node(), _kernel(kernel)
{
	_loopnode = 0;
}

void OuterLoopEnd::run(Stack<Manager::T_STACKITEM> & nodestack, bool usestack)
{
	if (_exitnode == 0 || _loopnode== 0) return;

	unsigned int size = _jobs.size();

	for (unsigned int i=0; i<size; i++)
	{
		unsigned int & stid = _jobs.data()[i];
#ifdef VERBOSE
		cout << "outerloopend " << stid << endl;
#endif	
		bool loop = _kernel.outerLoopEndModule(stid);

		unsigned int * nextjobs;
		if (loop)
		{
#ifdef VERBOSE
			cout << "again" << endl;
#endif	
			nextjobs = _loopnode->jobQueue().alloc(1);
		} else
		{
#ifdef VERBOSE
			cout << "exit outer loop" << endl;
#endif	
			nextjobs = _exitnode->jobQueue().alloc(1);
		}
		*nextjobs = stid;
	}

	_jobs.clear();

	if (usestack) {
		Manager::T_STACKITEM & loopitem = nodestack.push();
		loopitem.node = _loopnode;
		loopitem.usestack = true;

		Manager::T_STACKITEM & finishitem = nodestack.push();
		finishitem.node = _exitnode;
		finishitem.usestack = true;
	}



}


}}
