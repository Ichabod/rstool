#include "list.h"
#include <iomanip>
#include <iostream>

using namespace std;

namespace OptChol
{

List::List(Node*node)
{
    first=node;
    last=node;
    length=1;
    exist=true;
    tempexist=false;
    superfirst=0;
    superlast=0;
    superlength=0;
    hashindex=0;

}

List::List()
{
    first=0;
    last=0;
    exist=false;
    length=0;
    tempexist=false;
    superfirst=0;
    superlast=0;
    superlength=0;
    hashindex=0;
}

List::~List()
{

    Node*dummy;
    Node*dummy_2;
    dummy=first;
    dummy_2=dummy;
    if (length!=0)
    {
        if (length==1) delete first;
        else
        {

            for (int i=0;i<length-1;i++)
            {
                dummy=dummy->next;
                delete dummy_2;
                dummy_2=dummy;

            }
            delete dummy;
        }
    }
}


void List::add(int number)
{
    if(length==0){
   exist=true;
    first=new Node(number);
    last=first;
    length=1;


    }
    else{
    last->next=new Node(number);
    last->next->before=last;
    last=last->next;
    length++;
    }

}

void List::addfirst(int number)
{
    exist=true;
    first=new Node(number);
    last=first;
    length=1;
}

int List::getlength()
{

    return length;

}

void List::setlength()
{

    length--;
}
 void List::setlength(int number){
     length=number;

 }

void List::eliminate(int number)
{

    Node*dummy;
    bool test;
    test=true;
    if (length!=0)
    {

        if (length==1)
        {
           delete first;
            first=0;
        }

        else
        {
            if (first->value==number)
            {
                dummy=first;
                first=first->next;
                delete dummy;
                dummy=0;
            }
            else
            {
                dummy=first;

                for (int i=1;i<length-1&&test==true;i++)
                {
                    dummy=dummy->next;
                    if (dummy->value==number)
                    {

                        dummy->next->before=dummy->before;
                        dummy->before->next=dummy->next;
                        delete dummy;
                        dummy=0;
                        test=false;
                    }

                }

                if (test==true)
                {
                    dummy=dummy->next;
                    if (dummy->value==number)
                    {
                        last=last->before;

                        delete dummy;
                        dummy=0;
                        test=false;
                    }
                    else cout <<"Fehler"<<endl;
                }

            }
        }
        length--;
    }
}
void List::eliminatesup(int number)
{

    Node*dummy;
    bool test;
    test=true;
    if (length!=0)
    {

        if (length==1)
        {
           delete first;
            first=0;
        }

        else
        {
            if (first->value==number)
            {
                dummy=first;
                first=first->next;
                delete dummy;
                dummy=0;
            }
            else
            {
                dummy=first;

                for (int i=1;i<length-1&&test==true;i++)
                {
                    dummy=dummy->next;
                    if (dummy->value==number)
                    {

                        dummy->next->before=dummy->before;
                        dummy->before->next=dummy->next;
                        delete dummy;
                        dummy=0;
                        test=false;
                    }
                }

                if (test==true)
                {
                    dummy=dummy->next;
                    if (dummy->value==number)
                    {
                        last=last->before;
                        delete dummy;
                        dummy=0;
                        test=false;
                    }
                    else cout <<"Fehler"<<endl;
                }
            }
        }
        length--;
    }
}


void List::ausgabe()
{
    Node*walk;
    walk=first;
    cout <<"lenght: "<<length<<endl<<endl;
    if (length!=0)
    {
        for (int i=0;i<length-1;i++)
        {
            cout <<walk->value+1<<endl<<endl;
            walk=walk->next;
        }
        cout <<walk->value+1<<endl<<endl;
    }
    cout <<endl<<endl;

}

void List::eleminate(Node*node)
{

    if(length==1){
        first=0;
        last=0;
    }

    else if (node==first)
    {
        first=first->next;
        first->before=0;
    }
    else if (node==last)
    {
        last=last->before;
        last->next=0;

    }
    else
    {
        node->before->next=node->next;
        node->next->before=node->before;
    }
    length--;
}

void List::eleminatelast(Node*node)
{
    if (length==1)
    {
        first=0;
        last=0;
    }
    else
    {
        last=last->before;
        last->next=0;
    }
    length--;

}

void List::sort(){
   int wert;
   Node*dummy;
   Node*dummy_2;
    if(length>1){
        dummy=first->next;
        dummy_2=dummy;
     for (int i = 2; i<= length ;i++){
     wert=dummy->value;
     while( dummy!=first && dummy->before->value > wert){
          dummy->value = dummy->before->value;
           dummy=dummy->before;

     }
        dummy->value=wert;
        dummy_2=dummy_2->next;
        dummy=dummy_2;

     }

    }
}

 void List::addsupernode(int number,List*list){

     if(superlength==0) {
         superfirst=superlast=new Node(number);


     }
     else{
         superlast->next=new Node(number);
         superlast=superlast->next;
         }

         superlength++;

     if(list->superlength!=0){
         if(list->superlength==1) {
             superlast->next=list->superfirst;
             superlast=superlast->next;
             }

         else{
             superlast->next=list->superfirst;
             superlast=list->superlast;
         }
         superlength=superlength+list->superlength;
     }
   this->reallength=this->reallength+1+list->superlength;
     list->superlength=0;
     list->superfirst=0;
     list->superlast=0;
 }

    
}
