#ifndef _SCHEDULEJOBQUEUE_H
#define _SCHEDULEJOBQUEUE_H

#include <cstring>
#include <iostream>
#include <cstdlib>

namespace Scheduler
{
	template <typename T>
	class JobQueue
	{
		public:
			JobQueue()
			{
				_size=1;
				_pos=0;
				_queue = new T[_size];
			}

			~JobQueue()
			{
				delete[] _queue;
			}

			inline T * alloc(unsigned int requiredsize)
			{
				if (_pos+requiredsize > _size)
				{
//					std::cout << "realloc:" << _size << " -> " << _pos+requiredsize << std::endl;
					_size = _pos+requiredsize;
					T * temp;

					temp = new T[_size];
					if (_pos>0) memcpy(temp, _queue, sizeof(T)*_pos);

					delete[] _queue;

					_queue = temp;
				}
				unsigned int oldpos = _pos;
				_pos += requiredsize;

				return _queue+oldpos;
			}

			inline void discard(unsigned int size)
			{
				if (_pos > size)
				{
					_pos -= size;
				} else
				{
					_pos = 0;
				}
			}

			inline void importJobs(T * data, unsigned int size)
			{
				T * ptr = alloc(size);
				memcpy(ptr, data, sizeof(T)*size);
			}

			inline void exportJobs(T * data)
			{
				memcpy(data, _queue, sizeof(T)*_pos);
			}

			inline void clear()
			{
				_pos=0;
			}

			inline unsigned int size() const 
			{
				return _pos;
			}

			inline T * data()
			{
				return _queue;
			}

		protected:
			T * _queue;
			unsigned int _size;
			unsigned int _pos;

	};

}

#endif