#ifndef DATA_STORAGE_H
#define DATA_STORAGE_H

//#define STORAGE_SIZE 10
#define STORAGE_SIZE 10

#include <iostream>

namespace Scheduler
{
	class Storage
	{
	public:
		Storage()
		{
			for (unsigned int i=0; i<STORAGE_SIZE; i++)
			{
				_freecells[i] = i;	
			}
			
			_cursor = STORAGE_SIZE;
		}

		inline unsigned int alloc()
		{
			if (_cursor == 0) return 0xFFFFFFFF;

			--_cursor;

			return _freecells[_cursor];
		}

		inline unsigned int * alloc(unsigned int size)
		{
			if (_cursor < size) return 0;

			_cursor -= size;

			return &_freecells[_cursor];
		}

		inline void free(unsigned int cell)
		{
			if (_cursor < STORAGE_SIZE)	
			{
				_freecells[_cursor] = cell;
				++_cursor;
			}
		}		

		inline unsigned int countFreeCells() const
		{
			return _cursor;
		}

		inline unsigned int size() const
		{
			return STORAGE_SIZE-countFreeCells();
		}

		void printcontent()
		{
			for (unsigned int i=0; i<_cursor; i++)
			{
				std::cout << _freecells[i] << " ";
			}
			std::cout << std::endl;
		}


	protected:
		unsigned int _freecells[STORAGE_SIZE];
		unsigned int _cursor;

	};


}


#endif
