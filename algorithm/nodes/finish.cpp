#include "finish.h"
#include "../initialguess.h"
#include <iostream>
#include <cmath>

using namespace Scheduler;
using namespace std;

//#define VERBOSE

namespace Algorithm { namespace Nodes
{


Finish::Finish(Algorithm::Kernel & kernel, Algorithm::Buffer & buffer, Discr::Grid & grid, Scheduler::Storage & storage,
				pthread_cond_t & buffercond, pthread_mutex_t & condmutex, const int & threadnumber) :
	Node(), _kernel(kernel), _buffer(buffer), _grid(grid), _storage(storage), _buffercond(buffercond), _condmutex(condmutex), _threadnumber(threadnumber)
{
	devcounter=0;
	cellcounter=0;
}

void Finish::run(Stack<Manager::T_STACKITEM> & nodestack, bool usestack)
{
	unsigned int size = _jobs.size();

	_buffer.lock();
	int totaladded=0;

	for (unsigned int i=0; i<size; i++)
	{
		unsigned int & stid = _jobs.data()[i];

#ifdef VERBOSE
		cout << "ID " << stid << ": err= "<< _kernel.returncode[stid];
			cout << " outersteps=" << _kernel.outersteps[stid];
			cout << " steps=" << _kernel.steps[stid] << endl;
#endif

		_grid.setBusy(_kernel.grididx[stid], false);

		if (_kernel.returncode[stid]==Kernel::SUCCESS && !_grid.isMarked(_kernel.grididx[stid]))
		{		


			_grid.setMarked(_kernel.grididx[stid]);

//			if (devcounter < 10) // 
			{
				
				InitialGuess * guess = 0;
				int numcells=0;
				for (int j=0; j<_grid.countNeighbors(); j++)
				{
					int cellidx = _grid.getNeighborIdx(_kernel.grididx[stid], j);
					if (cellidx != -1 && !_grid.isMarked(cellidx))
					{
						numcells++;
						if (guess == 0)
						{
							double mu = max(min(_kernel.mu[stid]*1e4, 0.9), 1e-5);
							if (_kernel.Bcomputed[stid])
							{
								guess = new InitialGuess(0, _kernel.xi[stid], _kernel.vec_x[stid], _kernel.vec_z[stid],  _kernel.vec_y[stid],
									_kernel.s_b[stid], _kernel.y_b[stid], mu);	
							} else
							{
								guess = new InitialGuess(0, _kernel.xi[stid], _kernel.vec_x[stid], _kernel.vec_z[stid],  _kernel.vec_y[stid],
									0, 0, mu);	
							}
						}

						Buffer::ITEM & item = _buffer.push();
						item.guess = guess;
						item.grididx = cellidx;
					}
				}
				if (numcells > 0)
				{
					guess->count = numcells;
					totaladded+=numcells;
	#ifdef VERBOSE
					cout << "added " << numcells << " items" << endl;
	#endif
				}
			}


		}


		cellcounter++;
		_storage.free(stid);
	}

	_buffer.unlock();

	if (totaladded>0)
	{
	#ifdef VERBOSE
		cout << "[" << _threadnumber << "] added " << totaladded << ", broadcasting..." << endl;
	#endif

		if (_nextnode && usestack)
		{
			Manager::T_STACKITEM & startitem = nodestack.push();
			startitem.node = _nextnode;
			startitem.usestack = false;
		}
		pthread_mutex_lock(&_condmutex);
		pthread_cond_broadcast(&_buffercond);
		pthread_mutex_unlock(&_condmutex);
	}


	_jobs.clear();

	if (size > 0) devcounter++;

}


}}
