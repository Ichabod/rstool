
#ifndef ODEDOPRI853_H
#define ODEDOPRI853_H

#include "dopri.h"

namespace ODE
{

	class DoPri853 : public DoPri
	{

		public:
			DoPri853 ( ODEFUNC func, int ndgl );
	};

}

#endif
