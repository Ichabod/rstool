#ifndef OPTCHOLVARIABLES_H
#define OPTCHOLVARIABLES_H

/// Hilfklasse f�r Liste der aktiven Knoten

namespace OptChol
{
class Variables{

    public:
    Variables *next;
    bool first;
    Variables*before;
    int externaldegree;
    int index;
    Variables();
    ~Variables();

};

	
}

#endif