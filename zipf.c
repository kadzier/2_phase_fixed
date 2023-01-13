#include <math.h>
#include "zipf.h"
int dist_len; // the number of elements in the distribution
double udist[MAX_Z]; // underlying distribution
double sum_udist[MAX_Z]; // cumulative



void GenZipf(double alpha){
    double cursum = 0.0;
    for (int i=0;i<dist_len;i++){
        udist[i] = 1/pow(1.0+i, alpha);
        cursum += udist[i];
    }
    printf("%d\n",dist_len);
    for (int i=0;i<dist_len;i++){
        udist[i] /= cursum;
        sum_udist[i] = udist[i] + (i==0 ? 0 : sum_udist[i-1]);
        // printf("Dist[%d]: %lf, %lf\n", i, udist[i], sum_udist[i]);
    }
    return;
}
