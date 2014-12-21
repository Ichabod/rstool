#ifndef OPTCHOLQUOTIENTGRAPH_H
#define OPTCHOLQUOTIENTGRAPH_H

#include <set>
#include "list.h"
#include "variables.h"

namespace OptChol
{

class Quotientgraph{
    private:

    public:
     int dimension;
    Node **hashbucket;        /// Hashbucket
    List *A;                  /// Liste mit allen Adjazentlisten mit aktiven Knoten
    List *E;                   /// Liste mit allen Adjazenzlisten mit eliminierten Knoten
    Variables*variables;       /// Liste mit aktiven Knoten
    Variables*elements;        /// Liste mit allen eliminierten Knoten
    Variables*lastelement;     /// Zeiger auf den letzten eliminierten Knoten
    Variables*firstelement;    /// Zeiger auf ersten eliminierten Knoten
    Variables*firstvariable;   /// zeiger auf ersten noch aktiven Knoten

    Quotientgraph(int dimension, const std::set<unsigned int> & layout);
    ~Quotientgraph();

    inline static unsigned int coord_encode(unsigned int row, unsigned int col)
    {
        return (row << 16) + col;
    }

     void settrue(List* list);
     void setfalse(int number);
     void setfalseE(int number);
     void eliminationstep(int p);   ///Elimination von p aus dem Quotient-Graph
     void ausgabe();
     void crush(int number);
     void superremove(int p);                           /// Nicht unterscheidbare Knoten aus Graphen entfernen
     int hash(int number);
     bool indistinguistest(int number,int value);           /// testfunktion für nicht-unterscheidbarkeit
     int hashfind(int listindex);
     void resethash(int length);                             /// setzt Hashbucket auf NULL
     int findminimum();                                     /// Startknoten mit minimalem Grad
     void subset(int p);                                    /// Berechnung von L_e/L_p  für alle e
     void eliminatevariable(int index);                   /// aktive Knoten aus Liste entfernen
     void addelement(int index);                          ///  Eliminierte Knoten zur Liste hinzufügen
     void eliminateelement(int index);                     /// Eliminierte Knoten aus Liste entfernen
     int minexdegree();                                    /// Suche Knoten mit minimalem External Degree
     void exdegreeupdate(int p,int dimension,int counter);  /// Berechnung des AMD für alle relevanten Knoten
     void exdegreeausgabe();
};

    
}

#endif
