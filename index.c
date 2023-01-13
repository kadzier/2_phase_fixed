#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int TwoToOne(int x, int y){
    return (x+y)*(x+y+1)/2 + x;
}

long int LongTwoToOne(long int x, long int y){
    return (x+y)*(x+y+1)/2 + x;
}

void OneToTwo(int n, int *x, int *y){
    int xplusy = (int) floor(sqrt(2*n));  // [y(y+1)+2y]/2 = (y+2)(y+1);
   if (xplusy * (xplusy+1) > 2*n){
//        printf("* ");
        xplusy--;
    }
//   else printf("  ");
    *x = n - xplusy * (xplusy+1)/2;
    *y = xplusy - *x;
    return;
}


void LongOneToTwo(long int n, long int *x, long int *y){
    long int xplusy = (long int) floor(sqrt(2*n));  // [y(y+1)+2y]/2 = (y+2)(y+1);
   if (xplusy * (xplusy+1) > 2*n){
//        printf("* ");
        xplusy--;
    }
//   else printf("  ");
    *x = n - xplusy * (xplusy+1)/2;
    *y = xplusy - *x;
    return;
}


/* to fix the empty slots (due to no negative 0), if x or y is neg, increment
 by 1, store there instead (and still do the ++ or +2 offset).  On
 reversing, if negative value, negate and also subtract 1 */


int TwoToOneNeg(int x, int y){
    int idxbasic = TwoToOne(abs(x), abs(y)) * 4;
    if (x < 0){
        idxbasic++;
    }
    if (y < 0){
        idxbasic += 2;
    }
    return idxbasic;
}

void OneToTwoNeg(int n, int *x, int *y){
    OneToTwo(n/4,x,y);
    int negcheck = n - 0*4; // This needs fixing!!!
    if (negcheck > 1){
        *(y) *= -1;
        negcheck -= 2;
    }
    if (negcheck==1){
        *(x) *= -1;
    }
    return;
}
                            


int TwoToOneLimited(int x, int y){
    return x * (x+1) / 2 + y;
}

void OneToTwoLimited(int n, int *x, int *y){
    *x = (int) floor(sqrt(2*n));
    if (*x*(*x+1) > 2 *n){
        (*x)--;
    }
    *y = n - (*x * (*x+1) / 2);
    return;
}

long int LongTwoToOneLimited(long int x, long int y){
    return x * (x+1) / 2 + y;
}

void LongOneToTwoLimited(long int n, long int *x, long int *y){
    *x = (long int) floor(sqrt(2*n));
    if (*x*(*x+1) > 2 *n){
        (*x)--;
    }
    *y = n - (*x * (*x+1) / 2);
    return;
}

int TwoToOneLimited2(int x, int y){
    return TwoToOneLimited(x-1,y);
}

void OneToTwoLimited2(int n, int *x, int *y){
    OneToTwoLimited(n,x,y);
    (*x)++;
    return;
}

int NToOne(int N, long int *vals){
    long int curidx = vals[N-1];
    while (N>1){
        N--;
//        printf("(%d,%d) -> is ", vals[N-1], curidx);
        curidx = LongTwoToOne(vals[N-1], curidx);
//        printf("%d\n",curidx);
    }
    return curidx;
}

void OneToN(long int idx, int N, long int *valstore){
    int i=0;
    long int next_idx;
    while (i < N-1){
        LongOneToTwoLimited(idx, valstore+i, &next_idx);
//        printf("%d -> (%d,%d)\n", idx, valstore[i], next_idx);
        idx = next_idx;
        i++;
    }
    valstore[N-1] = idx;
    return;
}

int NToOneLimited(int N, long int *vals){
    long int curidx = vals[0];
    int i = 1;
    while (i<N){
        printf("(%ld,%ld) -> is ", vals[i], curidx);
        curidx = LongTwoToOneLimited(curidx, vals[i]);
        printf("%ld\n",curidx);
        i++;
    }
    return curidx;
}

void OneToNLimited (long int idx, int N, long int *valstore){
    int i=N-1;
    long int next_idx;
    while (i > 0){
        LongOneToTwoLimited(idx, &next_idx, valstore+i);
        printf("%ld -> (%ld,%ld)\n", idx, valstore[i], next_idx);
        idx = next_idx;
        i--;
    }
    valstore[0] = idx;
    return;
}

long int Choose(long int n, long int k){
    long int retval = 1;
    for (long int i=1;i<=k;i++){
        retval *= (n-k+i);
        retval /= i;
    }
    return retval;
}


long int ChooseList[100000];
int _choose_up_to = -1;


// special case: k <= n/2
int ChooseIdx(int n, int k){
    if (k*2>n){
        k=n-k;
    }
    return (n+1)*(n+1)/4 + k;
}

long int RChoose(int n, int k){
    int n_cur, k_cur, diff_cur;
    if (n<k){
        return 0;
    }
    if (2*k > n){
        return RChoose(n,n-k);
    }

    int idx = ChooseIdx(n,k);
    //printf("Storing at %d\n", idx);
    
    if (idx <= _choose_up_to){
        //printf("C(%d,%d)=%ld\n",n,k,ChooseList[idx]);
        return ChooseList[idx];
    }
    for (int i=_choose_up_to + 1; i <= idx; i++){
        n_cur =  (int)floor(2*sqrt((double)i)) - 1;
        k_cur = i - (n_cur+1)*(n_cur+1)/4;
        if (k_cur*2>n_cur){
            k_cur=0;
            n_cur++;
        }
//        printf("Doing idx %d = C(%d,%d)", i, n_cur, k_cur);
        if (k_cur==0){
            ChooseList[i] = 1;
        }
        else {
            ChooseList[i] = ChooseList[ChooseIdx(n_cur-1,k_cur)] +
                ChooseList[ChooseIdx(n_cur-1,k_cur-1)];
        }
//        printf("=%ld\n", ChooseList[i]);
    }
    _choose_up_to = idx;
    //printf("C(%d,%d)=%ld\n",n,k,ChooseList[idx]);
    
    return ChooseList[idx];
}

long int DecreasingTuple(int terms, long int *vals){
    long int retval = 0;
    for (int i=0;i<terms;i++){
        long int cval = Choose(vals[i], (long int) terms - i);
        //printf("cval=%ld\n", cval);
        retval += cval;
    }
    return retval;
}

void OneToDecreasing(long int the_one, int terms, long int *vals){
    for (int i=0;i<terms;i++){
        vals[i] = terms - i-1;
    }
    int cur_idx = 0;
    if (DecreasingTuple(terms, vals) == the_one){
        return;
    }
    while (cur_idx < terms){
        vals[cur_idx]++;
        printf("Trying "); for (int i=0;i<terms;i++){printf("%ld ", vals[i]);}printf("\n");
        if (DecreasingTuple(terms, vals) == the_one){
            return;
        }
        if (DecreasingTuple(terms,vals) >  the_one){
            vals[cur_idx]--;
            cur_idx++;
        }
    }
    printf("ERROR: could not do OneToDecreasing\n");
    exit(0);
}


long int NonIncreasingTuple(int terms, long int *vals){
    long int retval = 0;
    for (int i=0;i<terms;i++){
        long int cval = RChoose(vals[i]+terms-i-1, (long int) terms - i);
        retval += cval;
    }
    return retval;
}



void OneToNonIncreasing(long int the_one, int terms, long int *vals){
    for (int i=0;i<terms;i++){
        vals[i] = 0;
    }
    int cur_idx = 0;
    if (NonIncreasingTuple(terms, vals) == the_one){
        return;
    }
    while (cur_idx < terms){
        vals[cur_idx]++;
//        printf("Trying "); for (int i=0;i<terms;i++){printf("%ld ", vals[i]);}printf("\n");
        if (NonIncreasingTuple(terms, vals) == the_one){
            return;
        }
        if (NonIncreasingTuple(terms,vals) >  the_one){
            vals[cur_idx]--;
            cur_idx++;
        }
    }
    printf("ERROR: could not do OneToNonIncreasing\n");
    exit(0);
}


    






/*
int main(int argc, char **argv){
    int x = 1;
    int y = 0;
    int cursum=0;
    int xret, yret;
    int n;

    while (x<1000){
        n = TwoToOneLimited2(x,y);
        OneToTwoLimited2(n,&xret, &yret);
        printf("(%d,%d) -> %d -> (%d,%d)\n", x,y,n, xret,yret);

        x++;
        y--;
        if (y<0){
            y=x;
            x=0;
            printf("--------------------\n");
        }


        y++;
        if (y==x){
            y=0;
            x++;
        }
    }
    return 1;
}
        

*/

