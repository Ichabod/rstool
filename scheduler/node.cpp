#include "node.h"

#include "stack.h"

namespace Scheduler
{


Node::Node()
{
	_suggestedsize=1;
	_nextnode=0;
	_exitnode=0;
	_stopwatchsum=0.0;
}

}
