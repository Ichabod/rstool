#ifndef _SCHEDULERSTACK_H
#define _SCHEDULERSTACK_H

#include <cstring>

namespace Scheduler
{
	template <typename T>
	class Stack
	{
	public:
		Stack()
		{
			_memsize=1;
			_size=0;
			_stack = new T[_memsize];
		}

		~Stack()
		{
			delete[] _stack;
		}

		T & push()
		{
			if (_size == _memsize)
			{
				_memsize <<= 1;
				T * tmp = new T[_memsize];
				std::memcpy(tmp, _stack, sizeof(T)*_size);
				delete[] _stack;
				_stack = tmp;
			}
			return _stack[_size++];
		}
		T & pop()
		{
			return _stack[--_size];
		}
		unsigned int size()
		{
			return _size;
		}

	protected:
		T * _stack;
		unsigned int _size;
		unsigned int _memsize;


	};


}

#endif