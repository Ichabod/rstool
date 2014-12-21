#include "buffer.h"

#include <cstring>
#include <iostream>
#include <cstdlib>
#include <pthread.h>

using namespace std;

namespace Algorithm
{
   
Buffer::Buffer(int initialsize)
{
	_start = 0;
	_size = 0;
	_currentsize = initialsize;
	_buffer = new ITEM[_currentsize];

	pthread_mutex_init ( &_locker, NULL );
}

Buffer::~ Buffer()
{
   delete[] _buffer;
	pthread_mutex_destroy(&_locker);
}

void Buffer::printcontent()
{
   int __end = _start+_size;
   if (__end >= _currentsize)
   {
      __end -= _currentsize;
   }

   for (int i=0; i<_currentsize; i++)
   {
      if (i < __end  && __end <= _start && _size!=0) cout << _buffer[i].grididx;
      else if (i >= _start  && __end <= _start && _size!=0) cout << _buffer[i].grididx;
      else if (i >= _start && i < __end && __end > _start && _size!=0) cout << _buffer[i].grididx;
      else cout << ".";
   }

   cout << endl;
}

Buffer::ITEM & Buffer::push()
{
   int __end = _start+_size;
   if (__end >= _currentsize)
   {
      __end -= _currentsize;
   }

	if (_size == _currentsize)
	{
	   int destidx=_currentsize;
	   bool copy=true;
	   if (_start<__end)
	   {
	      copy=false;
	   }

      int oldsize = _currentsize;

      _currentsize *= 2;

      ITEM * tmpbuf = new ITEM[_currentsize];
      memcpy(tmpbuf, _buffer, sizeof(ITEM)*oldsize);
      delete[] _buffer;
      _buffer = tmpbuf;

      if (copy && __end != 0) memcpy(_buffer+destidx, _buffer, __end*sizeof(ITEM));
      __end=_start+_size;
	}

   _size ++;
   return _buffer[__end];
}

void Buffer::push(ITEM * src, int size)
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

      int oldsize = _currentsize;

      while (_size + size > _currentsize)
      {
         _currentsize *= 2;
      }
      ITEM * tmpbuf = new ITEM[_currentsize];
      memcpy(tmpbuf, _buffer, sizeof(ITEM)*oldsize);
      delete[] _buffer;
      _buffer = tmpbuf;

      if (copy && __end != 0) memcpy(_buffer+destidx, _buffer, __end*sizeof(ITEM));
      __end=_start+_size;
   }

   if (__end+size <_currentsize)
   {
      memcpy(&_buffer[__end], src, size*sizeof(ITEM));
   } else
   {
      memcpy(&_buffer[__end], src, (_currentsize - __end)*sizeof(ITEM));
      memcpy(&_buffer[0], &((char*)src)[_currentsize - __end], (size - (_currentsize-__end))*sizeof(ITEM));
   }
   _size += size;
}


Buffer::ITEM & Buffer::peek()
{
   return _buffer[_start];
}


int Buffer::popToBuffer(Buffer & dest, int size)
{
	if (size > _size)
	{
	   size=_size;
	}
	if (size==0) return 0;

   int __start = _start;
   if (__start+size <=_currentsize)
   {
      dest.push(_buffer+__start, size);
   } else
   {
      dest.push(_buffer+__start, _currentsize - __start);
      dest.push(_buffer, size - (_currentsize - __start));
   }

   discard(size);

   return size;
}


bool Buffer::discard(int size)
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

int Buffer::size()
{
	return _size;
}

void Buffer::clear()
{
	_size = 0;
	_start = 0;
}

void Buffer::lock()
{
	pthread_mutex_lock ( &_locker );
}

void Buffer::unlock()
{
	pthread_mutex_unlock ( &_locker );
}

}

