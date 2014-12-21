#include "innerloopend.h"

#include <iostream>
#include <cmath>

using namespace std;
using namespace Scheduler;

// #define VERBOSE

namespace Algorithm { namespace Nodes
{


InnerLoopEnd::InnerLoopEnd(Algorithm::Kernel & kernel) :
	Node(), _kernel(kernel)
{
	_loopnode = 0;
}

void InnerLoopEnd::run(Stack<Manager::T_STACKITEM> & nodestack, bool usestack)
{
	if (_nextnode == 0 || _loopnode== 0 || _exitnode == 0) return;

	unsigned int size = _jobs.size();

	for (unsigned int i=0; i<size; i++)
	{
		unsigned int & stid = _jobs.data()[i];
#ifdef VERBOSE
		cout << "innerloopend " << stid << endl;
#endif	
		Kernel::T_ERR err;
		bool loop = _kernel.innerLoopEndModule(stid, err);

		unsigned int * nextjobs;
		if (err != Kernel::SUCCESS)
		{
#ifdef VERBOSE
			cout << "error: " << err << endl;
#endif	
			nextjobs = _exitnode->jobQueue().alloc(1);
		} else if (loop)
		{
#ifdef VERBOSE
			cout << "again" << endl;
#endif	
			nextjobs = _loopnode->jobQueue().alloc(1);
		} else
		{
#ifdef VERBOSE
			cout << "exit inner loop" << endl;
#endif	
			nextjobs = _nextnode->jobQueue().alloc(1);
		}
		*nextjobs = stid;
	}

	if (usestack) {
		Manager::T_STACKITEM & nextitem = nodestack.push();
		nextitem.node = _nextnode;
		nextitem.usestack = true;
	}

	_jobs.clear();

}


}}
