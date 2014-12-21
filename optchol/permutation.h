#ifndef OPTCHOLPERMUTATION_H
#define OPTCHOLPERMUTATION_H

#include <set>

namespace OptChol
{

class Permutation
{

private:
    Permutation() {}

public:
    static void calc(int dimension, const std::set<unsigned int> & layout, int * permute);
    virtual ~Permutation() = 0;

};

    
}

#endif
