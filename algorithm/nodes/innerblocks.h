#ifndef _ALGORITHMNODESINNERBLOCKS_H
#define _ALGORITHMNODESINNERBLOCKS_H

#include "../../scheduler/node.h"
#include "../kernel.h"

namespace Algorithm { namespace Nodes
{
	class InnerBlocks : public Scheduler::Node
	{
		public:
			InnerBlocks(Algorithm::Kernel & kernel);

			void run(Scheduler::Stack<Scheduler::Manager::T_STACKITEM> & nodestack, bool usestack=true);
			
		protected:
			Algorithm::Kernel & _kernel;
	};


}}

#endif