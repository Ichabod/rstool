#ifndef _SCHEDULERNODE_H
#define _SCHEDULERNODE_H

#include "../tools/rtclock.h"
#include "jobqueue.h"
#include "stack.h"
#include "manager.h"

namespace Scheduler
{
	class Node
	{
		public:
			Node();

			virtual void run(Stack<Manager::T_STACKITEM> & nodestack, bool usestack=true) = 0;

			inline bool full() const
			{
				return _jobs.size() >= _suggestedsize;
			}

			inline unsigned int suggestedQueueSize() const
			{
				return _suggestedsize;
			}

			inline JobQueue<unsigned int> & jobQueue()
			{
				return _jobs;
			}
	
			inline void setNextNode(Node * nextnode)
			{
				_nextnode=nextnode;
			}

			inline void setExitNode(Node * exitnode)
			{
				_exitnode=exitnode;
			}

			inline void startStopWatch()
			{
				_stopwatch.reset();
			}

			inline void stopStopWatch()
			{
				_stopwatchsum += _stopwatch.elapsedSeconds();
			}

			inline double elapsedSeconds()
			{
				return _stopwatchsum;
			}



		protected:

			Node * _nextnode;
			Node * _exitnode;

			unsigned int _suggestedsize;
			JobQueue<unsigned int> _jobs;

			Tools::RTClock _stopwatch;
			double _stopwatchsum;

	};

}

#endif