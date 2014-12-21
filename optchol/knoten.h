#ifndef OPTCHOLKNOTEN_H
#define OPTCHOLKNOTEN_H

namespace OptChol
{


class Knoten{

    public:
    Knoten*next;
    Knoten*top;
    Knoten*down;
    Knoten();
    Knoten(int number,int value);
    int getcolumn();
     int valueindex;
    ~Knoten();
    int x;   ///zeilenindex
    int y;   ///spaltenindex
    private:









};
    
}

#endif