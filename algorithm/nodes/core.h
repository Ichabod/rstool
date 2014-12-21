#ifndef _ALGORITHMNODESCORE_H
#define _ALGORITHMNODESCORE_H

#include "../../scheduler/node.h"
#include "../kernel.h"

namespace Algorithm { namespace Nodes
{
	class Core : public Scheduler::Node
	{
		public:
			Core(Algorithm::Kernel & kernel);

			void run(Scheduler::Stack<Scheduler::Manager::T_STACKITEM> & nodestack, bool usestack=true);
			
		protected:
			Algorithm::Kernel & _kernel;
	};


}}

#endif