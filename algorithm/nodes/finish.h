#ifndef _ALGORITHMFINISHED_H
#define _ALGORITHMFINISHED_H

#include "../../scheduler/node.h"
#include "../../scheduler/storage.h"
#include "../buffer.h"
#include "../../discr/grid.h"
#include "../kernel.h"
#include <fstream>
#include <pthread.h>

namespace Algorithm { namespace Nodes
{

	class Finish : public Scheduler::Node
	{
		public:
			Finish(Algorithm::Kernel & kernel, Algorithm::Buffer & buffer, Discr::Grid & grid, Scheduler::Storage & storage,
				pthread_cond_t & buffercond, pthread_mutex_t & condmutex, const int & threadnumber);

			void run(Scheduler::Stack<Scheduler::Manager::T_STACKITEM> & nodestack, bool usestack=true);

			inline int getCellCounter() const
			{
				return cellcounter;
			}

		protected:
			Algorithm::Kernel & _kernel;
			Algorithm::Buffer & _buffer;
			Discr::Grid & _grid;
			Scheduler::Storage & _storage;
			pthread_cond_t & _buffercond;
			pthread_mutex_t & _condmutex;
			const int & _threadnumber;

			int devcounter;
			int cellcounter;

	};

}}

#endif