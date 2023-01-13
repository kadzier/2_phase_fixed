#ifndef _ZIPF_H
#define _ZIPF_H

#define MAX_Z 100000000
#define TL TwoToOneLimited

extern int dist_len; // the number of elements in the distribution
extern double udist[MAX_Z]; // underlying distribution
extern double sum_udist[MAX_Z]; // cumulative


void GenZipf(double alpha);
// run after GenZipf so hae Zipf dist which will majorize




#endif
