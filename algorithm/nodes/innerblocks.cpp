#include "innerblocks.h"
#include <stdio.h>
#include <iostream>

using namespace Scheduler;
using namespace std;

// #define VERBOSE

namespace Algorithm { namespace Nodes
{

InnerBlocks::InnerBlocks(Algorithm::Kernel & kernel) :
	Node(), _kernel(kernel)
{
}

void InnerBlocks::run(Stack<Manager::T_STACKITEM> & nodestack, bool usestack)
{
	if (_nextnode==0 || _exitnode==0) return;

	for (int j=0; j<_kernel.N; j++)
	{

		for (unsigned int i=0; i<_jobs.size(); i++)
		{
			unsigned int & stid = _jobs.data()[i];

			if (_kernel.returncode[stid] == Kernel::SUCCESS)
			{
				_kernel.returncode[stid] = _kernel.innerBlockUModule(stid, j);
			}
		}

		for (unsigned int i=0; i<_jobs.size(); i++)
		{
			unsigned int & stid = _jobs.data()[i];

			if (_kernel.returncode[stid] == Kernel::SUCCESS)
			{
				_kernel.returncode[stid] = _kernel.innerBlockXModule(stid, j);
			}
		}

		for (unsigned int i=0; i<_jobs.size(); i++)
		{
			unsigned int & stid = _jobs.data()[i];

			if (_kernel.returncode[stid] == Kernel::SUCCESS)
			{
				_kernel.returncode[stid] = _kernel.innerBlockMergeModule(stid, j);
			}
		}

	}

	for (unsigned int i=0; i<_jobs.size(); i++)
	{

		unsigned int & stid = _jobs.data()[i];

		unsigned int * nextjobs;
		if (_kernel.returncode[stid] == Kernel::SUCCESS)
		{
			nextjobs = _nextnode->jobQueue().alloc(1);
		} else
		{
			nextjobs = _exitnode->jobQueue().alloc(1);
		}
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
