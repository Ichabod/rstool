#ifndef _ALGORITHMNODESPREPARESIGRHS_H
#define _ALGORITHMNODESPREPARESIGRHS_H

#include "../../scheduler/node.h"
#include "../kernel.h"

namespace Algorithm { namespace Nodes
{
	class PrepareSigRhs : public Scheduler::Node
	{
		public:
			PrepareSigRhs(Algorithm::Kernel & kernel);

			void run(Scheduler::Stack<Scheduler::Manager::T_STACKITEM> & nodestack, bool usestack=true);
			
		protected:
			Algorithm::Kernel & _kernel;
	};


}}

#endif