#include "rtclock.h"
#include <cstring>
#include <cmath>

namespace Tools
{


RTClock::RTClock()
{
	gettimeofday ( & ( _starttime ), NULL );
}

RTClock::~ RTClock()
{
}

double RTClock::elapsedSeconds()
{
	timeval end;
	gettimeofday ( &end, NULL );

	double result = ( double ) end.tv_sec- ( double ) _starttime.tv_sec + ( double ) ( end.tv_usec - _starttime.tv_usec )  /1000000.0;

	return result;

}

void RTClock::reset()
{
	gettimeofday ( & ( _starttime ), NULL );
}


void RTClock::reset ( const timeval & starttime )
{
	_starttime.tv_sec = starttime.tv_sec;
	_starttime.tv_usec = starttime.tv_usec;
}

void RTClock::getStartTime ( timeval & time )
{
	time.tv_sec = _starttime.tv_sec;
	time.tv_usec = _starttime.tv_usec;
}

void RTClock::setElapsedSeconds(double seconds)
{
	timeval time;
	gettimeofday ( &time, NULL );
	time.tv_sec -= (time_t)floor(seconds);
	time_t usec = (time_t)((seconds-floor(seconds))*1E+6);

	if (time.tv_usec < usec)
	{
		time.tv_usec += 1E+6;
		time.tv_sec -= 1;
	}
	time.tv_usec -= usec;

	_starttime.tv_sec = time.tv_sec;
	_starttime.tv_usec = time.tv_usec;

}

}
