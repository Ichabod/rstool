#include "mergeblocks.h"
#include <stdio.h>
#include <iostream>

using namespace Scheduler;
using namespace std;

// #define VERBOSE

namespace Algorithm { namespace Nodes
{

MergeBlocks::MergeBlocks(Algorithm::Kernel & kernel) :
	Node(), _kernel(kernel)
{
}

void MergeBlocks::run(Stack<Manager::T_STACKITEM> & nodestack, bool usestack)
{
	if (_nextnode==0 || _exitnode==0) return;

	for (unsigned int i=0; i<_jobs.size(); i++)
	{

		unsigned int & stid = _jobs.data()[i];
#ifdef VERBOSE
		cout << "MergeBlocks " << stid << endl;
#endif	
		Kernel::T_ERR err = _kernel.mergeBlocksModule(stid);

		unsigned int * nextjobs;
		if (err == Kernel::SUCCESS)
		{
			nextjobs = _nextnode->jobQueue().alloc(1);
		} else
		{
			nextjobs = _exitnode->jobQueue().alloc(1);
		}
		_kernel.returncode[stid] = err;
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
