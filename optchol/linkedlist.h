#ifndef OPTCHOLLINKEDLIST_H
#define OPTCHOLLINKEDLIST_H

#include "knoten.h"
#include "list.h"

#include <set>

namespace OptChol
{

 /// Listenstruktur für L-Matrix

class Linkedlist{

    public:
    Knoten*first;  ///Zeiger auf ersten Knoten
    Knoten*last;  /// Zeiger auf letzten Knoten
    Linkedlist();
    ~Linkedlist();

    void add(int number,int value);
    Knoten* add(List*list,int counter,List*A,Knoten *end);
    void ausgabe();
     void ausgabeL();
     int ausgabe_2();
    void addcolumn(int number,int value,Knoten *end);
    void sort();
    void copy(int number,Linkedlist*list,int index);
    void ausgabe_3();
};
    
}

#endif