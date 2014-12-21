#ifndef ODEDOPRI5_H
#define ODEDOPRI5_H

#include "dopri.h"

namespace ODE
{


	class DoPri5 : public DoPri
	{

		public:
			DoPri5 ( ODEFUNC func, int ndgl );

	};

}

#endif
