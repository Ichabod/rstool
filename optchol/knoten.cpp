#include "knoten.h"


using namespace std;

namespace OptChol
{

Knoten::Knoten(){

    next=0;
    top=0;
    down=0;
    x=0;
    y=0;
    valueindex=-1;
}

Knoten::Knoten(int number,int value){

    next=0;
    top=0;
    down=0;
    x=number;
    y=value;
    valueindex=-1;

}

Knoten::~Knoten(){


}

int Knoten::getcolumn(){

    return y;

}
    
}
