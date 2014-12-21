#include "node.h"


namespace OptChol
{

Node::Node(int number){

    value=number;
     next=0;
    before=0;
      hashnext=0;
}
 Node:: Node(){


    value=0;
     next=0;
    before=0;
    hashnext=0;
 }
Node::~Node(){



}
	
}
