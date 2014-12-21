#ifndef OPTCHOLLIST_H
#define OPTCHOLLIST_H

#include "node.h"

namespace OptChol
{

/// Listenstruktur f�r A_i und E_i im Quotienten Graphen
/// Jedes listenobjekt A_i enth�lt auch alle Informationen �ber den Knoten i


class List{

 private:
        int length;   /// Anzahl der Nachbarknoten


 public:

    int memoryindex;   ///Hilsvariable f�r Hashfunktion
    int hashindex;    ///Hilsvariable f�r Hashfunktion
    int reallength;   /// Exakter Grad
    Node*superfirst;  /// Zeiger auf ersten nicht unterscheidbaren Knoten der mit i zusammengefasst wurde
    Node*superlast;
    int superlength;  /// Anzahl der anliegenden Nicht-unterscheidbaren Knoten
    Node*first;        /// Zeiger auf ersten Knoten in der Liste
    Node*last;
    bool tempexist;    /// Hilfvariable f�r Mischen von Listen
    bool exist;         /// Hilfvarible f�r Listenoperationen
    List();
    List(Node* node);
    ~List();
    void add(int number);
    int getlength();
    void setlength();
    void setlength(int number);
    void addfirst(int number);                     /// Ersten Knoten hinzuf�gen
    void eliminate(int number);                    ///Knoten loeschen �ber index
    void ausgabe();
    void eleminate(Node*node);                   ///Knoten l�schen �ber Zeiger
    void eleminatelast(Node*node);                   /// letzten Knoten l�schen
    void addsupernode(int number,List*list);      /// Nicht unterschiedbare Knoten zusammenfassen
    void sort();
    void eliminatesup(int number);                /// Nicht unterscheibare Knoten entfernen
};
    
}

#endif