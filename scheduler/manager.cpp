#include "manager.h"
#include "stack.h"
#include "node.h"

namespace Scheduler
{

Manager::Manager(Node & startnode) : _startnode(startnode)
{

}

void Manager::exec()
{
	Stack<T_STACKITEM> nodestack;

	T_STACKITEM & startitem = nodestack.push();
	startitem.node = &_startnode;
	startitem.usestack = true;

	while ( nodestack.size() > 0 )
	{
		T_STACKITEM currentnode = nodestack.pop();

		currentnode.node->startStopWatch();
		currentnode.node->run(nodestack, currentnode.usestack);
		currentnode.node->stopStopWatch();
	}
}


}