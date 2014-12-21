#ifndef TOOLSFIFOBUFFER_H
#define TOOLSFIFOBUFFER_H

#include <pthread.h>

namespace Tools
{

	class FIFOBuffer
		{
			public:
				FIFOBuffer(int initialsize);
				~FIFOBuffer();

				void push(void * src, int size);
				int peek(void * dest, int size);
				int peek(void * dest, int position, int size);
				int pop(void * dest, int size);
				int discard(int size);
				int popToFIFO(FIFOBuffer * buffer, int size);

				int size();

				void clear();

				void lock();
				void unlock();

				void printcontent();

			protected:

				inline int end(int size);

				pthread_mutex_t _locker;

				char * _buffer;
				int _start;
				int _size;
				int _currentsize;
		};
}
#endif // FIFOBUFFER_H
