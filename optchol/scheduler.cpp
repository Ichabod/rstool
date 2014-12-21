#include "scheduler.h"

#include <list>
#include <set>
#include <iostream>
#include <cmath>
#include <cstring>
#include "quotientgraph.h"
#include "permutation.h"

using namespace std;

// #define VERBOSE

namespace OptChol
{

int * Scheduler::generate(int n, double * A)
{
	int counter;

	int * packedinfo;

	int * n_rhselim = new int[n-1];
	int * n_matrixelim = new int[n-1];
	int * n_backsubst = new int[n-1];
	int * sdiag = new int[n];
	list<int> rhselim;
	list<int> backsubst;
	list<int> matrixelim;

	set<unsigned int> nonzeros;

	for (int j=0; j<n; j++)
	{
		for (int i=0; i<n; i++)
		{
			if (abs(*A) != 0.0)
			{
				nonzeros.insert(OptChol::Quotientgraph::coord_encode(i,j));
#ifdef VERBOSE				
				cout << "1 ";
			} else
			{
				cout << "0 ";
#endif
			}
			A++;
		}

#ifdef VERBOSE				
		cout << endl;
#endif
	}

	OptChol::Permutation::calc(n, nonzeros, sdiag);
#ifdef VERBOSE				
	cout << "PERMUTATION: ";
#endif
	for (int i=0; i<n; i++)
	{
		sdiag[i]++; // Fortran idx
		sdiag[i] = i+1;
#ifdef VERBOSE				
		cout << sdiag[i] << " ";
#endif
	}
#ifdef VERBOSE				
	cout << endl;
#endif

	int globalcounter=0;

	for (int i=0; i<n;i++)
	{
#ifdef VERBOSE				
		cout << "Dec Schritt " << i << endl;
#endif

		int i2 = sdiag[i]-1;
   		nonzeros.insert(OptChol::Quotientgraph::coord_encode(i2,i2));
#ifdef VERBOSE				
		cout << "  g(" << i2 << "," << i2 << ")" << endl;
#endif
      
		if (i < n-1)
		{
			counter = 0;
			for (int j=i+1; j<n; j++)
			{
				int j2 = sdiag[j]-1;
				int i2 = sdiag[i]-1;
				if (nonzeros.find(OptChol::Quotientgraph::coord_encode(j2,i2)) != nonzeros.end())
				{
#ifdef VERBOSE				
					cout << "  b(" << j2 << ") -= b(" << i2 << ")*g(" << j2 << "," << i2 << ")" << endl;
#endif
					rhselim.push_back( j2 );
					counter++;
				}

			}
			n_rhselim[i] = counter;

			counter = 0;
			for (int j=i+1; j <n; j++)
			{
				for (int k=i+1; k<=j; k++)
				{
					int i2 = sdiag[i]-1;
					int j2 = sdiag[j]-1;
					int k2 = sdiag[k]-1;
					if (nonzeros.find(OptChol::Quotientgraph::coord_encode(j2,i2)) != nonzeros.end() && nonzeros.find(OptChol::Quotientgraph::coord_encode(k2,i2)) != nonzeros.end())
					{
						if (nonzeros.find(OptChol::Quotientgraph::coord_encode(j2,k2)) == nonzeros.end())
						{
							nonzeros.insert(OptChol::Quotientgraph::coord_encode(j2,k2));
#ifdef VERBOSE				
							cout << "fillin: g(" << j2 << "," << k2 << ")" << endl;
#endif
						}
#ifdef VERBOSE				
						cout << "  g(" << j2 << "," << k2 << ") -= g(" << j2 << "," << i2 << ")*g(" << k2 << "," << i2 << ")" << endl;
#endif
						matrixelim.push_back( j2 );
						matrixelim.push_back( k2 );
						counter++;
					}
					globalcounter++;
				}
			}
			n_matrixelim[i] = counter;
      }
   }

   for (int i=0; i<n-1; i++)
   {
#ifdef VERBOSE				
		cout << "BS Schritt " << i << endl;
#endif

		counter = 0;
		int bidx = n-i-1;
		for (int j=0; j<bidx; j++)
		{
			int bidx2 = sdiag[bidx]-1;
			int j2 = sdiag[j]-1;
			if (nonzeros.find(OptChol::Quotientgraph::coord_encode(bidx2,j2)) != nonzeros.end())
			{
#ifdef VERBOSE				
				cout << "  b(" << j2 << ") -= b(" << bidx2 << ")*g(" << bidx2 << "," << j2 << ")" << endl;
#endif
				backsubst.push_back( j2 );
				counter++;
			}
		}

		n_backsubst[i] = counter;
   	}
#ifdef VERBOSE				
   	cout << "Alle eliminationen: " << globalcounter << endl;
   	cout << "Sparse: " << matrixelim.size()/2 << endl;
#endif

	int totalsize = 7 + n + (n-1) + (n-1) + (n-1) + rhselim.size() + matrixelim.size() + backsubst.size();
	packedinfo = new int[totalsize];
	packedinfo[0] = 8;
	packedinfo[1] = packedinfo[0]+n;
	packedinfo[2] = packedinfo[1]+(n-1);
	packedinfo[3] = packedinfo[2]+(n-1);
	packedinfo[4] = packedinfo[3]+(n-1);
	packedinfo[5] = packedinfo[4]+rhselim.size();
	packedinfo[6] = packedinfo[5]+matrixelim.size();
	memcpy(packedinfo+packedinfo[0]-1, sdiag, sizeof(int)*n);
	memcpy(packedinfo+packedinfo[1]-1, n_rhselim, sizeof(int)*(n-1));
	memcpy(packedinfo+packedinfo[2]-1, n_matrixelim, sizeof(int)*(n-1));
	memcpy(packedinfo+packedinfo[3]-1, n_backsubst, sizeof(int)*(n-1));

	int * ptr = packedinfo+packedinfo[4]-1;
	for (list<int>::iterator it = rhselim.begin(); it != rhselim.end(); ++it)
	{
		*ptr = (*it)+1;
		ptr++;
	}

	for (list<int>::iterator it = matrixelim.begin(); it != matrixelim.end(); ++it)
	{
		*ptr = (*it)+1;
		ptr++;
	}

	for (list<int>::iterator it = backsubst.begin(); it != backsubst.end(); ++it)
	{
		*ptr = (*it)+1;
		ptr++;
	}

	delete[] n_rhselim;
	delete[] n_matrixelim;
	delete[] n_backsubst;
	delete[] sdiag;

	return packedinfo;
}



}