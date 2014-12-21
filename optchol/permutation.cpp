#include "permutation.h"
#include "linkedlist.h"
#include "quotientgraph.h"

using namespace std;

namespace OptChol
{

void Permutation::calc(int dimension, const std::set<unsigned int> & layout, int * permute)
{
    int length;

    Quotientgraph quotientgraph(dimension, layout);

    Linkedlist*L=new Linkedlist[dimension];
    Node*dummy;
    Node*walk;
    int steps;
    int index=0;
    Knoten*L_dummy;
    int counter=0;
    int counterdummy=0;
    for (int i=0;i<dimension;i++) L[i].add(i,i);

    index=quotientgraph.findminimum();
    while (counter<dimension)
    {
        permute[counter]=index;
        quotientgraph.addelement(index);

        quotientgraph.eliminationstep(index);

        length= quotientgraph.hashfind(index);

        quotientgraph.eliminatevariable(index);
        quotientgraph.resethash(length);

        L_dummy=L[counter].first;

        /// Erzeugen der Spalte in L
        L_dummy=L[counter].add(&(quotientgraph.A[index]),counter,quotientgraph.A,L_dummy);

        steps=quotientgraph.A[index].superlength;
        counterdummy=counter;
        counter++;



/// Nicht-unterscheidbare Knoten eliminieren

        if (steps!=0)
        {
            walk=quotientgraph.A[index].superfirst;
            quotientgraph.A[index].superfirst=0;
            dummy=walk;

            for (int i=0;i<steps-1;i++)
            {
                L[counterdummy].addcolumn(walk->value,counterdummy,L_dummy);
                L_dummy=L_dummy->down;
                walk=walk->next;
            }
            L[counterdummy].addcolumn(walk->value,counterdummy,L_dummy);

            walk=dummy;

            for (int i=0;i<steps-1;i++)
            {
                L[counter].copy(walk->value,&L[counter-1],counter);
                permute[counter]=walk->value;
                walk=walk->next;
                delete dummy;
                dummy=walk;
                counter++;
            }

            L[counter].copy(walk->value,&L[counter-1],counter);
            permute[counter]=walk->value;
            delete dummy;
            counter++;
        }

        if (counter<dimension)
        {
            quotientgraph.exdegreeupdate(index,dimension,counter);
            index=quotientgraph.minexdegree();
        }

        //quotientgraph.exdegreeausgabe();
    }

    delete[] L;
}

}




