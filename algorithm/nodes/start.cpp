#include "start.h"
#include "../fortranheader.h"

//#define VERBOSE

using namespace Scheduler;

namespace Algorithm { namespace Nodes
{

Start::Start(Algorithm::Kernel & kernel, Algorithm::Buffer & buffer, Discr::Grid & grid, Scheduler::Storage & storage,
			bool * & threaddone, pthread_cond_t & buffercond, pthread_mutex_t & condmutex, const int & threadnumber, const int & threadcount, int & statstoragesize) :
	Node(), _kernel(kernel), _buffer(buffer), _grid(grid), _storage(storage),
	_threaddone(threaddone), _buffercond(buffercond), _condmutex(condmutex), _threadnumber(threadnumber), _threadcount(threadcount), _statstoragesize(statstoragesize)
{
}

Start::~Start()
{
}

void Start::run(Stack<Manager::T_STACKITEM> & nodestack, bool usestack)
{
	bool abort = false;

	_buffer.lock();
	if (_storage.size() == 0 && _buffer.size() == 0)
	{
#ifdef VERBOSE
		printf("[%d] empty...\n", _threadnumber);
#endif
		_buffer.unlock();

		abort=true;

		pthread_mutex_lock(&_condmutex);
		_threaddone[_threadnumber] = true;
		for (int i=0; i<_threadcount; i++)
		{
			if (!_threaddone[i])
			{
				abort=false;
				break;
			}
		}

		if (!abort)
		{
#ifdef VERBOSE
			printf("[%d] wait, some other threads are still running...\n", _threadnumber);
#endif
			pthread_cond_wait(&_buffercond, &_condmutex);
			pthread_mutex_unlock(&_condmutex);

			_buffer.lock();
			if (_buffer.size()>0)
			{
				_threaddone[_threadnumber]=false;
			}
		} else
		{
#ifdef VERBOSE
			printf("[%d] all threads done, goodbye...\n", _threadnumber);
#endif
			pthread_cond_broadcast(&_buffercond);
			pthread_mutex_unlock(&_condmutex);

			_buffer.lock();
		}

	}

	int refillcounter = 0;
	int requeuecounter = 0;
	if (!abort && _buffer.size() > 0 && _nextnode != 0)
	{
		int prebuffersize = _buffer.size();
		int bufferdeccounter = prebuffersize;
		while (bufferdeccounter>prebuffersize-std::max(prebuffersize/_threadcount,1) && _storage.countFreeCells() > 0)
		{
			Buffer::ITEM & bufferitem = _buffer.peek();

			if (!_grid.isBusy(bufferitem.grididx))
			{
				if (!_grid.isMarked(bufferitem.grididx))
				{
					_grid.setBusy(bufferitem.grididx, true);
					unsigned int * job = _storage.alloc(1);
					_nextnode->jobQueue().importJobs(job, 1);

					unsigned int & stid = *job;
					
					_kernel.grididx[stid]=bufferitem.grididx;
					_grid.readCellBounds(bufferitem.grididx, _kernel.lb_x[stid], _kernel.ub_x[stid]);

					_kernel.initModule(bufferitem.guess, stid);
					refillcounter++;
				}
				InitialGuess::Destroy(bufferitem.guess);
			} else
			{
				Buffer::ITEM & requeueitem = _buffer.push();
				requeueitem = bufferitem;
				requeuecounter++;
			}
			_buffer.discard(1);	
			bufferdeccounter--;		
		}

#ifdef VERBOSE
			printf("[%d] refill %d items...\n", _threadnumber, refillcounter);
			if (requeuecounter>0)
				printf("[%d] requeue %d busy items...\n", _threadnumber, requeuecounter);
#endif
	}

	if (usestack && !abort)
	{
		Manager::T_STACKITEM & thisitem = nodestack.push();
		thisitem.node = this;
		thisitem.usestack=true;
		if (refillcounter>0)
		{
			Manager::T_STACKITEM & nextitem = nodestack.push();
			nextitem.node = _nextnode;
			nextitem.usestack=true;
		}
	}

	_statstoragesize = _storage.size();
	_buffer.unlock();
/*
	if (refillcounter == 0 && !abort)
	{
#ifdef VERBOSE
		printf("[%d] wait, no new cells refilled...\n", _threadnumber);
#endif
		pthread_mutex_lock(&_condmutex);
		pthread_cond_wait(&_buffercond, &_condmutex);
		pthread_mutex_unlock(&_condmutex);
	}
*/

}

}}