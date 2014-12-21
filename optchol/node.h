#ifndef OPTCHOLNODE_H
#define OPTCHOLNODE_H

///  Knotenklasse für Listen A_i

namespace OptChol
{

class Node{

    private:

    public:
    Node *next;
    Node*hashnext;  ///Zeiger auf nächstes Knoten im Hashbucket
    Node*before;
    int value;
    Node(int number);
    ~Node();
    Node();
};

	
}

#endif