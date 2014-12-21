#include "fifobuffer.h"

#include <cstring>
#include <iostream>
#include <cstdlib>
#include <pthread.h>

using namespace std;

namespace Tools
{
   
FIFOBuffer::FIFOBuffer(int initialsize)
{
	_start = 0;
	_size = 0;
	_currentsize = initialsize;
	_buffer = (char*)malloc(_currentsize);

	pthread_mutex_init ( &_locker, NULL );
}

FIFOBuffer::~ FIFOBuffer()
{
	free(_buffer);
	pthread_mutex_destroy(&_locker);
}

void FIFOBuffer::printcontent()
{
   int __end = _start+_size;
   if (__end >= _currentsize)
   {
      __end -= _currentsize;
   }

   for (int i=0; i<_currentsize; i++)
   {
      if (i < __end  && __end <= _start && _size!=0) cout << _buffer[i];
      else if (i >= _start  && __end <= _start && _size!=0) cout << _buffer[i];
      else if (i >= _start && i < __end && __end > _start && _size!=0) cout << _buffer[i];
      else cout << ".";
   }

   cout << endl;
}

void FIFOBuffer::push(void * src, int size)
{
   int __end = _start+_size;
   if (__end >= _currentsize)
   {
      __end -= _currentsize;
   }

	if (_size+size > _currentsize)
	{
	   int destidx=_currentsize;
	   bool copy=true;
	   if (_start<__end)
	   {
	      copy=false;
	   }

      while (_size + size > _currentsize)
      {
         _currentsize *= 2;
      }
      _buffer=(char*)realloc(_buffer, _currentsize);

      if (copy && __end != 0) memcpy(&_buffer[destidx], &_buffer[0], __end);
      __end=_start+_size;
	}

   if (__end+size <_currentsize)
   {
      memcpy(&_buffer[__end], src, size);
   } else
   {
      memcpy(&_buffer[__end], src, _currentsize - __end);
      memcpy(&_buffer[0], &((char*)src)[_currentsize - __end], size - (_currentsize-__end));
   }
   _size += size;
}

int FIFOBuffer::peek(void * dest, int position, int size)
{
	if (position+size > _size)
	{
	   size=_size-position;
	}
	if (size<0) return 0;

   int __start = _start+position;
   if (__start+size <=_currentsize)
   {
      memcpy(dest, &_buffer[__start], size);
   } else
   {
      memcpy(dest, &_buffer[__start], _currentsize - __start);
      memcpy(&((char*)dest)[_currentsize - __start], &_buffer[0], size - (_currentsize - __start));
   }

   return size;
}

int FIFOBuffer::popToFIFO(FIFOBuffer * dest, int size)
{
	if (size > _size)
	{
	   size=_size;
	}
	if (size==0) return 0;

   int __start = _start;
   if (__start+size <=_currentsize)
   {
      dest->push(&_buffer[__start], size);
   } else
   {
      dest->push(&_buffer[__start], _currentsize - __start);
      dest->push(&_buffer[0], size - (_currentsize - __start));
   }

   discard(size);

   return size;
}

int FIFOBuffer::peek(void * dest, int size)
{
	return peek(dest, 0, size);
}

int FIFOBuffer::discard(int size)
{
   if (size > _size) size=_size;

   _start += size;
   if (_start >= _currentsize)
   {
      _start -= _currentsize;
   }
   _size -= size;

   return size;
}

int FIFOBuffer::pop(void * dest, int size)
{
   int readbytes = peek(dest, size);
   if (readbytes>0) discard(size);

   return readbytes;
}

int FIFOBuffer::size()
{
	return _size;
}

void FIFOBuffer::clear()
{
	_size = 0;
	_start = 0;
}

int FIFOBuffer::end(int size)
{
	int __end = _start+_size+size;
	return __end<_currentsize ? __end : (__end - _currentsize);
}

void FIFOBuffer::lock()
{
	pthread_mutex_lock ( &_locker );
}

void FIFOBuffer::unlock()
{
	pthread_mutex_unlock ( &_locker );
}

}

