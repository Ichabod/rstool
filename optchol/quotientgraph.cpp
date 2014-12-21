#include "quotientgraph.h"

#include <iostream>

using namespace std;

namespace OptChol
{

Quotientgraph::Quotientgraph(int dimension, const std::set<unsigned int> & layout)
{
    this->dimension=dimension;
    hashbucket=new Node*[dimension+1];
    for (int i=0;i<dimension+1;i++) hashbucket[i]=0;
    A=new List[dimension]();
    E=new List[dimension]();
    variables=new Variables[dimension];
    elements=new Variables[dimension];
    firstelement=0;

    for (set<unsigned int>::iterator it = layout.begin(); it != layout.end(); it++)
    {
        int row = (*it) >> 16;
        int col = (*it) & 0x0000FFFF;
        if (row>col)
        {
            if (A[row].exist==true) A[row].add(col);
            else A[row].addfirst(col);

            if (A[col].exist==true) A[col].add(row);
            else A[col].addfirst(row);
        }
    }
    for (int i=0;i<dimension;i++) A[i].reallength=A[i].getlength();

    variables[0].first= true;
    variables[0].externaldegree=A[0].getlength();
    variables[0].index=0;
    if(dimension>1){
     variables[0].next=&variables[1];
    for (int i=1;i<dimension-1;i++)
    {
        variables[i].next=&variables[i+1];
        variables[i].before=&variables[i-1];
        variables[i].externaldegree=A[i].getlength();
        variables[i].index=i;
    }
    variables[dimension-1].before=&variables[dimension-2];
    variables[dimension-1].externaldegree=A[dimension-1].getlength();
    variables[dimension-1].next=0;
    variables[dimension-1].index=dimension-1;

    }
      firstvariable=&variables[0];
}
Quotientgraph::~Quotientgraph()
{
    delete[]A;
    delete[]E;
    delete[]hashbucket;
    delete[]variables;
    delete[]elements;

}
void Quotientgraph::settrue(List* list)
{
    Node*walk;
    int dummy;
    dummy=list->getlength();
    if (dummy!=0)
    {
        walk=list->first;
        for (int i=0;i<dummy-1;i++)
        {
            A[walk->value].tempexist=true;
            walk=walk->next;
        }
        A[walk->value].tempexist=true;
    }
}

void Quotientgraph::setfalse(int number)
{
    Node*walk;
    int dummy;
    dummy=A[number].getlength();
    if (dummy!=0)
    {
        walk=A[number].first;
        for (int i=0;i<dummy-1;i++)
        {
            A[walk->value].tempexist=false;
            walk=walk->next;
        }
        A[walk->value].tempexist=false;
    }
}
void Quotientgraph::setfalseE(int number)
{

    Node*walk;
    int dummy;
    dummy=E[number].getlength();
    if (dummy!=0)
    {
        walk=E[number].first;
        for (int i=0;i<dummy-1;i++)
        {
            A[walk->value].tempexist=false;
            walk=walk->next;
        }
        A[walk->value].tempexist=false;
    }
}

void Quotientgraph::eliminationstep(int p)
{
    A[p].tempexist=true;
    A[p].exist=false;
    Node*walk;
    Node*walk_2;
    int dummy;
    int dummy_2;
    int dummy_3;
    Node*help;
    help=0;

    settrue(&A[p]);

    ///L_P=A_P...
    if (E[p].getlength()!=0)
    {

        walk=E[p].first;

        for (int i=0;i<E[p].getlength()-1;i++)
        {
            dummy=walk->value;
            dummy_3=A[dummy].getlength();


            if ( dummy_3!=0)
            {
                walk_2=A[dummy].first;
                for (int j=0;j<dummy_3-1;j++)
                {
                    dummy_2= walk_2->value;

                    if (A[dummy_2].tempexist==false)
                    {
                        A[dummy_2].tempexist=true;

                        A[p].add(dummy_2);
                        A[p].reallength=A[p].reallength+1+A[dummy_2].superlength;
                    }

                    walk_2=walk_2->next;
                }
                dummy_2=walk_2->value;

                if (A[dummy_2].tempexist==false)
                {
                    A[dummy_2].tempexist=true;
                    A[p].add(dummy_2);
                    A[p].reallength=A[p].reallength+1+A[dummy_2].superlength;
                }
            }
            walk=walk->next;
        }

        dummy=walk->value;
        dummy_3=A[dummy].getlength();
        if (dummy_3!=0)
        {
            walk_2=A[dummy].first;

            for (int j=0;j<dummy_3-1;j++)
            {
                dummy_2= walk_2->value;
                if (A[dummy_2].tempexist==false)
                {
                    A[dummy_2].tempexist=true;
                    A[p].add(dummy_2);
                    A[p].reallength=A[p].reallength+1+A[dummy_2].superlength;

                }
                walk_2=walk_2->next;
            }
            dummy_2= walk_2->value;
            if (A[dummy_2].tempexist==false)
            {
                A[dummy_2].tempexist=true;
                A[p].add(dummy_2);
                A[p].reallength=A[p].reallength+1+A[dummy_2].superlength;
            }
        }
    }

    /// A_i=(A_i\L_p)\p

    dummy=A[p].getlength();
    if (dummy!=0)
    {
        walk=A[p].first;
        for (int i=0;i<dummy-1;i++)
        {
            dummy_2=walk->value;
            if (E[dummy_2].getlength()==0)  E[dummy_2].addfirst(p);
            else E[dummy_2].add(p);
            walk_2=A[dummy_2].first;
            dummy_3=A[dummy_2].getlength();

            if (dummy_3!=0)
            {
                for (int j=0;j<dummy_3-1;j++)
                {
                    if (A[walk_2->value].tempexist==true)
                    {
                        A[dummy_2].eleminate(walk_2);
                        A[dummy_2].reallength=A[dummy_2].reallength-1-A[walk_2->value].superlength;
                        help=walk_2;
                        walk_2=walk_2->next;
                        delete help;
                        help=0;
                    }
                    else  walk_2=walk_2->next;
                }
                if (A[walk_2->value].tempexist==true)
                {
                    A[dummy_2].eleminatelast(walk_2);
                    A[dummy_2].reallength=A[dummy_2].reallength-1-A[walk_2->value].superlength;
                    delete walk_2;
                    walk_2=0;
                }
            }
            walk=walk->next;
        }

        dummy_2=walk->value;
        if (E[dummy_2].getlength()==0)  E[dummy_2].addfirst(p);
        else E[dummy_2].add(p);
        dummy_3=A[dummy_2].getlength();

        if (dummy_3!=0)
        {
            walk_2=A[dummy_2].first;
            for (int j=0;j<dummy_3-1;j++)
            {
                if (A[walk_2->value].tempexist==true)
                {
                    A[dummy_2].eleminate(walk_2);
                    A[dummy_2].reallength=A[dummy_2].reallength-1-A[walk_2->value].superlength;
                    help=walk_2;
                    walk_2=walk_2->next;
                    delete help;
                    help=0;
                }
                else walk_2=walk_2->next;
            }
            if (A[walk_2->value].tempexist==true)
            {
                A[dummy_2].eleminatelast(walk_2);
                A[dummy_2].reallength=A[dummy_2].reallength-1-A[walk_2->value].superlength;
                delete walk_2;
                walk_2=0;
            }
        }
    }
    setfalse(p);
    A[p].tempexist=false;



    dummy=E[p].getlength();

    if (dummy!=0)
    {
        walk=E[p].first;
        for (int j=0;j<dummy-1;j++)
        {
            A[walk->value].eliminate(p);
            walk=walk->next;
        }
        A[walk->value].eliminate(p);
    }

    if (dummy!=0)
    {
        walk=E[p].first;
        walk_2=walk;
        for (int i=0;i<dummy-1;i++)
        {
            walk_2=walk->next;
            eliminateelement(walk->value);
            crush(walk->value);
            delete walk;
            walk=0;
            E[p].setlength();
            walk=walk_2;
        }

        eliminateelement(walk->value);
        crush(walk_2->value);
        delete walk_2;
        walk_2=0;
        E[p].setlength();
    }
    E[p].exist=false;
}

void Quotientgraph::ausgabe()
{
    for (int i=0;i<dimension;i++)
    {
        cout<<"A["<<i+1<<"] : ";
        A[i].ausgabe();
        cout<<endl<<endl;
        cout<<"E["<<i+1<<"] : ";

        if (E[i].getlength()!=0)
        {
            E[i].ausgabe();
        }
        cout<<endl<<endl;
    }
}

void Quotientgraph::crush(int number)
{
    Node*walk;
    Node*del;
    int dummy;
    dummy=A[number].getlength();

    if (dummy!=0)
    {
        walk=A[number].first;
        del =walk;

        for (int i=0;i<dummy-1;i++)
        {
            E[walk->value].eliminate(number);
            walk=walk->next;
            A[number].eleminate(del);
            delete del;
            del=walk;
        }
        E[walk->value].eliminate(number);
        A[number].eleminate(del);

        delete walk;
        walk=0;
    }
}

void Quotientgraph::superremove(int p)
{
    int dummy;
    Node*walk;
    Node*help;
    dummy=E[p].getlength();

    if (dummy!=0)
    {
        walk=E[p].first;

        for (int j=0;j<dummy-1;j++)
        {
            A[walk->value].eliminatesup(p);
            help=walk->next;
            delete walk;
            walk=help;
        }

        A[walk->value].eliminatesup(p);
        delete walk;

        E[p].setlength(0);
    }

    dummy=A[p].getlength();
    if (dummy!=0)
    {

        walk=A[p].first;
        for (int j=0;j<dummy-1;j++)
        {
            A[walk->value].eliminatesup(p);
            help=walk->next;
            delete walk;
            walk=help;
        }
        A[walk->value].eliminatesup(p);
        delete walk;
        A[p].setlength(0);
    }
}

int Quotientgraph::hash(int number)
{
    int count =0;
    Node*walk;
    int  dummy=A[number].getlength();

    if (dummy!=0)
    {
        walk=A[number].first;
        for (int i=0;i<dummy-1;i++)
        {
            count=1+count+walk->value;
            walk=walk->next;
        }
        count=1+count+walk->value;

    }
    dummy=E[number].getlength();

    if (dummy!=0)
    {
        walk=E[number].first;
        for (int i=0;i<dummy-1;i++)
        {
            count=1+count+walk->value;
            walk=walk->next;
        }
        count=1+count+walk->value;
    }
    count=count%(dimension-1);
    count++;

    return count;
}
bool Quotientgraph::indistinguistest(int number,int value)
{
    bool test=true;
    if (A[number].getlength()!=A[value].getlength()|| E[number].getlength()!=E[value].getlength()) test=false;
    else
    {
        Node *walk;
        this->settrue(&A[number]);
        int dummy=A[value].getlength();
        if (dummy!=0)
        {
            walk=A[value].first;

            for (int i=0;i<dummy-1&&test==true;i++)
            {
                if (A[walk->value].tempexist==false)  test=false ;

                walk=walk->next;
            }
            if (A[walk->value].tempexist==false) test=false;
        }
        this->setfalse(number);

        this->settrue(&E[number]);
        dummy=E[value].getlength();
        if (dummy!=0&&test==true)
        {
            walk=E[value].first;

            for (int i=0;i<dummy-1&&test==true;i++)
            {
                if (A[walk->value].tempexist==false)  test=false ;
                walk=walk->next;
            }
            if (A[walk->value].tempexist==false) test=false;
        }
        this->setfalseE(number);
    }

    return test;
}


int Quotientgraph::hashfind(int listindex)
{
    bool test;
    bool test_2;
    int hashcount=0;
    int dummy;
    Node*walk;
    Node*go;
    Node*check;

    int steps=A[listindex].getlength();
    int superremovecount=0;

    if (steps!=0)
    {
        walk=A[listindex].first;
        for (int i=0;i<steps-1;i++)
        {
            test_2=true;
            dummy=hash(walk->value);

            if (hashbucket[dummy]==0)
            {
                hashbucket[dummy]=walk;
                A[hashcount].hashindex=dummy;
                hashcount++;
                walk=walk->next;
            }

            else
            {
                check=go= hashbucket[dummy];
                while (go!=0&&test_2==true)
                {
                    test= this->indistinguistest(walk->value,go->value);

                    if (test==true)
                    {
                        A[superremovecount].memoryindex=walk->value;
                        superremovecount++;
                        A[go->value].addsupernode(walk->value,&A[walk->value]);
                        test_2=false;
                    }

                    else go=go->hashnext;
                    if (go!=0) check=go;
                }

                if (test_2==true) check->hashnext=walk;
                walk=walk->next;

            }
        }

        dummy=hash(walk->value);
        test_2=true;

        if (hashbucket[dummy]==0)
        {
            hashbucket[dummy]=walk;
            A[hashcount].hashindex=dummy;
            hashcount++;
        }

        else
        {
            check=go= hashbucket[dummy];
            while (go!=0&&test_2==true)
            {
                test=this-> indistinguistest(walk->value,go->value);
                if (test==true)
                {
                    A[superremovecount].memoryindex=walk->value;
                    superremovecount++;
                    A[go->value].addsupernode(walk->value,&A[walk->value]);
                    test_2=false;
                }
                else go=go->hashnext;
                if (go!=0) check=go;
            }

            if (test_2==true) check->hashnext=walk;
        }

        for (int i=0;i<superremovecount;i++)
        {
            this->superremove(A[i].memoryindex);
            eliminatevariable(A[i].memoryindex);
        }
    }
    return hashcount;
}

void Quotientgraph::resethash(int length)
{
    Node*walk;
    Node*dummy;
    for (int i=0;i<length;i++)
    {
        walk=hashbucket[A[i].hashindex];
        hashbucket[A[i].hashindex]=0;

        while (walk!=0)
        {
            dummy=walk->hashnext;
            walk->hashnext=0;
            walk=dummy;
        }
    }
}
int Quotientgraph::findminimum()
{
    int dummy=0;
    int index=0;
    int min=A[0].getlength();

    for (int i=1;i<dimension;i++)
    {
        dummy=A[i].getlength();
        if (dummy<min)
        {
            min =dummy;
            index=i;
        }
    }
    return index;

}
int Quotientgraph::minexdegree()
{
    int index;
    int dummy=0;
    Variables*variable=firstvariable;
    index=variable->index;
    dummy=firstvariable->externaldegree;
    while (variable->next!=0)
    {
        variable=variable->next;
        if (variable->externaldegree<dummy)
        {
            dummy=variable->externaldegree;
            index=variable->index;
        }
    }
    return index;
}
void Quotientgraph::subset(int p)
{

    int degree_2;
    Node*walk;
    Node*walk_2;

    Variables*walk_3=firstelement;
    while (walk_3!=0)
    {
        A[walk_3->index].superlength=-1;
        walk_3=walk_3->next;
    }

    int degree=A[p].getlength();
    if (degree!=0)
    {
        walk=A[p].first;

        for (int i=0;i<degree-1;i++)
        {
            degree_2=E[walk->value].getlength();
            if (degree_2!=0)
            {
                walk_2=E[walk->value].first;

                for (int j=0;j<degree_2-1;j++)
                {
                    if (A[walk_2->value].superlength==-1)
                    {
                        A[walk_2->value].superlength=A[walk_2->value].reallength;
                    }

                    A[walk_2->value].superlength=A[walk_2->value].superlength-A[walk->value].superlength-1;

                    walk_2=walk_2->next;
                }
                if (A[walk_2->value].superlength==-1)
                {
                    A[walk_2->value].superlength=A[walk_2->value].reallength;
                }

                A[walk_2->value].superlength=A[walk_2->value].superlength-A[walk->value].superlength-1;
            }

            walk=walk->next;
        }
        degree_2=E[walk->value].getlength();
        if (degree_2!=0)
        {
            walk_2=E[walk->value].first;

            for (int j=0;j<degree_2-1;j++)
            {
                if (A[walk_2->value].superlength==-1)
                {
                    A[walk_2->value].superlength=A[walk_2->value].reallength;
                }
                A[walk_2->value].superlength=A[walk_2->value].superlength-A[walk->value].superlength-1;
                walk_2=walk_2->next;
            }

            if (A[walk_2->value].superlength==-1)
            {
                A[walk_2->value].superlength=A[walk_2->value].reallength;
            }
            A[walk_2->value].superlength=A[walk_2->value].superlength-A[walk->value].superlength-1;
        }
    }
}
void Quotientgraph::exdegreeupdate(int p,int dimension,int counter)
{
    int exdegree=0;
    int dummy=0;
    int degree_2=0;
    Node*walk;
    Node*walk_2;
    subset(p);
    int degree=A[p].getlength();
    if (degree!=0)
    {
        walk=A[p].first;

        for (int i=0;i<degree-1;i++)
        {

            exdegree=dimension-counter;
            if (variables[walk->value].externaldegree+A[p].reallength-A[walk->value].superlength-1<exdegree)  exdegree=variables[walk->value].externaldegree+A[p].reallength-A[walk->value].superlength-1;

            dummy=A[walk->value].reallength-A[walk->value].superlength+A[p].reallength-A[walk->value].superlength-1;
            degree_2=E[walk->value].getlength();

            if (degree_2!=0)
            {
                walk_2=E[walk->value].first;
                for (int j=0;j<degree_2-1;j++)
                {
                    if (walk_2->value!=p)
                    {
                        if (A[walk_2->value].superlength<0) dummy=dummy+A[walk_2->value].reallength;
                        else dummy=dummy+ A[walk_2->value].superlength;
                    }
                    walk_2=walk_2->next;

                }
                if (walk_2->value!=p)
                {
                    if (A[walk_2->value].superlength<0) dummy=dummy+A[walk_2->value].reallength;
                    else dummy=dummy+ A[walk_2->value].superlength;
                }
            }

            if (dummy<exdegree)
            {
                exdegree=dummy;
            }
            variables[walk->value].externaldegree=exdegree;
            walk=walk->next;
        }
        exdegree=dimension-counter;

        if (variables[walk->value].externaldegree+A[p].reallength-A[walk->value].superlength-1<exdegree)  exdegree=variables[walk->value].externaldegree+A[p].reallength-A[walk->value].superlength-1;

        dummy=A[walk->value].reallength-A[walk->value].superlength+A[p].reallength-A[walk->value].superlength-1;

        degree_2=E[walk->value].getlength();
        if (degree_2!=0)
        {

            walk_2=E[walk->value].first;
            for (int j=0;j<degree_2-1;j++)
            {
                if (walk_2->value!=p)
                {
                    if (A[walk_2->value].superlength<0) dummy=dummy+A[walk_2->value].reallength;
                    else dummy=dummy+ A[walk_2->value].superlength;
                }
                walk_2=walk_2->next;
            }

            if (walk_2->value!=p)
            {
                if (A[walk_2->value].superlength<0) dummy=dummy+A[walk_2->value].reallength;
                else dummy=dummy+ A[walk_2->value].superlength;
            }
        }

        if (dummy<exdegree)
        {
            exdegree=dummy;
        }
        variables[walk->value].externaldegree=exdegree;
    }
}

void Quotientgraph::eliminatevariable(int index)
{

    if (variables[index].first==false)
    {
        if (variables[index].next!=0)
        {
            variables[index].next->before= variables[index].before;
            variables[index].before->next= variables[index].next;
        }
        else
        {
            variables[index].before->next=0;
        }
    }
    else
    {
        if ( variables[index].next!=0)
        {
            firstvariable=variables[index].next;
            variables[index].next->first= true;

        }

        else {}
    }
}
void Quotientgraph::eliminateelement(int index)
{
    if (&elements[index]==firstelement)
    {
        if (firstelement==lastelement) firstelement=lastelement=0;

        else
        {
            firstelement=firstelement->next;
            firstelement->before=0;
        }
    }
    else if (&elements[index]==lastelement)
    {
        lastelement=lastelement->before;
        lastelement->next=0;
    }

    else
    {
        elements[index].before->next= elements[index].next;
        elements[index].next->before= elements[index].before;
    }

}
void Quotientgraph::addelement(int index)
{
    if (firstelement==0)
    {
        elements[index].index=index;
        firstelement=lastelement=&elements[index];
    }

    else
    {
        elements[index].index=index;
        lastelement->next=&elements[index];
        lastelement->next->before=lastelement;
        lastelement=lastelement->next;
    }
}
void Quotientgraph::exdegreeausgabe()
{
    Variables*walk=firstvariable;
    while (walk!=0)
    {
        cout<<walk->externaldegree<<endl;
        walk=walk->next;
    }

    cout<<endl;
}


    
}
