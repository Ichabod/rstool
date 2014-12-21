#ifndef _ALGORITHMNODESTART_h
#define _ALGORITHMNODESTART_h

#include <pthread.h>

#include "../../scheduler/node.h"
#include "../../scheduler/storage.h"
#include "../buffer.h"
#include "../../discr/grid.h"
#include "../kernel.h"

namespace Algorithm { namespace Nodes {


	class Start : public Scheduler::Node
	{
	public:
		Start(Algorithm::Kernel & kernel, Algorithm::Buffer & buffer, Discr::Grid & grid, Scheduler::Storage & storage,
			bool * & threaddone, pthread_cond_t & buffercond, pthread_mutex_t & condmutex, const int & threadnumber, const int & threadcount, int & statstoragesize);
		~Start();

		void run(Scheduler::Stack<Scheduler::Manager::T_STACKITEM> & nodestack, bool usestack=true);

	protected:
		Algorithm::Kernel & _kernel;
		Algorithm::Buffer & _buffer;
		Discr::Grid & _grid;
		Scheduler::Storage & _storage;
		bool * &_threaddone;
		pthread_cond_t & _buffercond;
		pthread_mutex_t & _condmutex;
		const int & _threadnumber;
		const int & _threadcount;
		int & _statstoragesize;

	};



}}


#endif