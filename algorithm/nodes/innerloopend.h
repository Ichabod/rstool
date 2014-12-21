#ifndef _ALGORITHMNODESINNERLOOPEND_H
#define _ALGORITHMNODESINNERLOOPEND_H

#include "../../scheduler/node.h"
#include "../kernel.h"
#include <fstream>

namespace Algorithm { namespace Nodes
{

	class InnerLoopEnd : public Scheduler::Node
	{
		public:
			InnerLoopEnd(Algorithm::Kernel & kernel);

			void run(Scheduler::Stack<Scheduler::Manager::T_STACKITEM> & nodestack, bool usestack=true);

			void setLoopNode(Node * loopnode)
			{
				_loopnode = loopnode;
			}

		protected:
			Algorithm::Kernel & _kernel;
			
			Scheduler::Node * _loopnode;

	};

}}

#endif