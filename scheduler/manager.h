#ifndef _SCHEDULERMANAGER_H
#define _SCHEDULERMANAGER_H

namespace Scheduler
{
	class Node;
	
	class Manager
	{
	public:
		typedef struct
		{
			Node * node;
			bool usestack;
		} T_STACKITEM;

		Manager(Node & startnode);

		void exec();


	protected:
		Node & _startnode;

	};


}

#endif