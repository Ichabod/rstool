#include "linkedlist.h"

#include <iostream>
#include <set>
using namespace std;

namespace OptChol
{


Linkedlist::Linkedlist()
{
    first=0;
    last=0;
}

Linkedlist::~Linkedlist()
{
    Knoten*dummy;
    while (first!=0)
    {
        dummy=first;

        first=first->down;
        delete dummy;
        dummy=first;
    }
}

void Linkedlist::add(int number,int value)
{

    if (first==0)
    {
        first=new Knoten(number,value);
        last=first;
        last->next=0;

    }
    else
    {
        Knoten*dummy=first;
        first=new Knoten(number,value);
        first->next=dummy;
    }

}
void Linkedlist::ausgabe()
{

    Knoten*dummy=first;

    while (dummy!=last)
    {

        cout<<dummy->valueindex<<" ";
        dummy=dummy->next;

    }
    cout<<dummy->valueindex<<" ";
    cout <<endl<<endl;
}
void Linkedlist::ausgabeL()
{

    Knoten*dummy=first;

    while (dummy!=last)
    {

        cout<<dummy->getcolumn()<<" ";
        dummy=dummy->next;
    }
    cout<<dummy->getcolumn()<<" ";

    cout <<endl<<endl;
}
int Linkedlist::ausgabe_2()
{
    int counter=0;
    Knoten*dummy=first;

    while (dummy!=0)
    {

        counter++;
        dummy=dummy->down;

    }
    return counter;
    cout <<endl<<endl;
}
void Linkedlist::ausgabe_3()
{
    Knoten*dummy=first;

    while (dummy!=0)
    {
        cout <<dummy->x<<endl;
        dummy=dummy->down;
    }
}
Knoten* Linkedlist::add(List*list,int counter,List*A,Knoten *end)
{

    int steps;
    int steps_2;
    Node*walk;
    Node*walk_2;
    Knoten*least;
    least=end;

    steps=list->getlength();

    if (steps!=0)
    {
        walk=list->first;
        for (int i=0;i<steps-1;i++)
        {

            this->addcolumn(walk->value,counter,least);
            least=least->down;
            steps_2=A[walk->value].superlength;

            if (steps_2!=0)
            {
                walk_2=A[walk->value].superfirst;
                for (int j=0;j<steps_2-1;j++)
                {

                    this->addcolumn(walk_2->value,counter,least);
                    least=least->down;
                    walk_2=walk_2->next;
                }
                this->addcolumn(walk_2->value,counter,least);
                least=least->down;
            }
            walk=walk->next;
        }
        this->addcolumn(walk->value,counter,least);

        least=least->down;
        steps_2=A[walk->value].superlength;

        if (steps_2!=0)
        {
            walk_2=A[walk->value].superfirst;
            for (int j=0;j<steps_2-1;j++)
            {
                this->addcolumn(walk_2->value,counter,least);
                least=least->down;
                walk_2=walk_2->next;
            }
            this->addcolumn(walk_2->value,counter,least);
            least=least->down;
        }
    }
    return least;
}

void Linkedlist::addcolumn(int number,int value,Knoten *end)
{
    Knoten*dummy=end;
    end->down=new Knoten(number,value);
    end->down->top=dummy;

}
void Linkedlist::sort()
{

    int wert;
    Knoten*walk=first->down;
    Knoten*walk_2;
    Knoten*one;
    one =walk;

    if (walk!=0)
    {
        walk=walk->down;

        while (walk!=0)
        {
            walk_2=walk;
            wert=walk_2->x;
            while (walk_2!=one&&wert<walk_2->top->x)
            {
                walk_2->x=walk_2->top->x;
                walk_2=walk_2->top;
            }
            walk_2->x=wert;
            walk=walk->down;
        }
    }
}
void Linkedlist::copy(int number,Linkedlist*list,int index)
{

    Knoten*walk=this->first;
    Knoten*walk_2=list->first;
    walk_2=walk_2->down;
    while (walk_2!=0)
    {

        if (walk_2->x!=number)
        {
            walk->down=new Knoten(walk_2->x,index);
            walk->down->top=walk;
            walk=walk->down;
        }
        walk_2=walk_2->down;
    }
}

}