#ifndef _ALGORITHMNODESMERGEBLOCKS_H
#define _ALGORITHMNODESMERGEBLOCKS_H

#include "../../scheduler/node.h"
#include "../kernel.h"

namespace Algorithm { namespace Nodes
{
	class MergeBlocks : public Scheduler::Node
	{
		public:
			MergeBlocks(Algorithm::Kernel & kernel);

			void run(Scheduler::Stack<Scheduler::Manager::T_STACKITEM> & nodestack, bool usestack=true);
			
		protected:
			Algorithm::Kernel & _kernel;
	};


}}

#endif