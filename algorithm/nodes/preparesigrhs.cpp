#include "preparesigrhs.h"
#include <stdio.h>
#include <iostream>

using namespace Scheduler;
using namespace std;

// #define VERBOSE

namespace Algorithm { namespace Nodes
{

PrepareSigRhs::PrepareSigRhs(Algorithm::Kernel & kernel) :
	Node(), _kernel(kernel)
{
}

void PrepareSigRhs::run(Stack<Manager::T_STACKITEM> & nodestack, bool usestack)
{
	if (_nextnode==0) return;

	for (unsigned int i=0; i<_jobs.size(); i++)
	{

		unsigned int & stid = _jobs.data()[i];
#ifdef VERBOSE
		cout << "PrepareSigRhs " << stid << endl;
#endif	
		_kernel.prepareSigRhsModule(stid);

		unsigned int * nextjobs;
		nextjobs = _nextnode->jobQueue().alloc(1);
		*nextjobs = stid;
	}

	if (usestack && _jobs.size() > 0) {
		Manager::T_STACKITEM & nextitem = nodestack.push();
		nextitem.node = _nextnode;
		nextitem.usestack = true;
	}

	_jobs.clear();

}

}}
