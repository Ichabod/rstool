#ifndef OPTCHOLFORTRANHEADER_H
#define OPTCHOLFORTRANHEADER_H

extern "C"
{
	void chlsco_(int * n, int * m, int * sb, int * rb, double * A, double * B, int * err, int * pinfo);
	void rrslvo_(int * n, int * m, int * sb, int * rb, double * A, double * B, int * pinfo);
}


#endif