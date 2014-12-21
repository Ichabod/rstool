#ifndef UTILSRTCLOCK_H
#define UTILSRTCLOCK_H

#include <sys/time.h>

namespace Tools
{
	class RTClock
	{

		public:
			RTClock();

			virtual ~RTClock();

			void reset();
			void reset( const timeval & starttime );

			double elapsedSeconds();
			void setElapsedSeconds(double seconds);
			void getStartTime ( timeval & time );


		protected:
			timeval _starttime;		/**< @brief Initial time for timer */

	};


}

#endif
