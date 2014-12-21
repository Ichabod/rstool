#ifndef ALGORITHMBUFFER_H
#define ALGORITHMBUFFER_H

#include <pthread.h>

namespace Algorithm
{
	class InitialGuess;

	class Buffer
	{
		public:
			typedef struct
			{
				InitialGuess * guess;
				int grididx;	
			} ITEM;

			Buffer(int initialsize);
			~Buffer();

			ITEM & push();
			void push(ITEM * src, int size);
			ITEM & peek();
			bool discard(int size);
			int popToBuffer(Buffer & buffer, int size);

			int size();

			void clear();

			void lock();
			void unlock();

			void printcontent();

		protected:

			pthread_mutex_t _locker;

			ITEM * _buffer;
			int _start;
			int _size;
			int _currentsize;
	};
}
#endif // ALGORITHMBUFFER_H
