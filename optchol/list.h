#ifndef OPTCHOLLIST_H
#define OPTCHOLLIST_H

#include "node.h"

namespace OptChol
{

/// Listenstruktur für A_i und E_i im Quotienten Graphen
/// Jedes listenobjekt A_i enthält auch alle Informationen über den Knoten i


class List{

 private:
        int length;   /// Anzahl der Nachbarknoten


 public:

    int memoryindex;   ///Hilsvariable für Hashfunktion
    int hashindex;    ///Hilsvariable für Hashfunktion
    int reallength;   /// Exakter Grad
    Node*superfirst;  /// Zeiger auf ersten nicht unterscheidbaren Knoten der mit i zusammengefasst wurde
    Node*superlast;
    int superlength;  /// Anzahl der anliegenden Nicht-unterscheidbaren Knoten
    Node*first;        /// Zeiger auf ersten Knoten in der Liste
    Node*last;
    bool tempexist;    /// Hilfvariable für Mischen von Listen
    bool exist;         /// Hilfvarible für Listenoperationen
    List();
    List(Node* node);
    ~List();
    void add(int number);
    int getlength();
    void setlength();
    void setlength(int number);
    void addfirst(int number);                     /// Ersten Knoten hinzufügen
    void eliminate(int number);                    ///Knoten loeschen über index
    void ausgabe();
    void eleminate(Node*node);                   ///Knoten löschen über Zeiger
    void eleminatelast(Node*node);                   /// letzten Knoten löschen
    void addsupernode(int number,List*list);      /// Nicht unterschiedbare Knoten zusammenfassen
    void sort();
    void eliminatesup(int number);                /// Nicht unterscheibare Knoten entfernen
};
    
}

#endif