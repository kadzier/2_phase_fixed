#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "index.c"
#include "zipf.c"
#include "index.c"
#include <sys/time.h>

#define _NEWVER

#define MAXBLOOMSIZE 500000000
#define EPSILON 0.00001
#define ETA_START_JUMP 0.001

// to terminate early if will take too long
#define MAX_HOUR_EST_RUN 15
#define MIN_SEC_MAKE_CALL 60
#define MIN_JL_MAKE_CALL 1000000
#define EST_ALPHA 0.9
#define DELTA_MAX 0.1
#define MAX_K 1000

FILE *thefile;

typedef struct {
    int BloomSize;
    long long Sigma;
    long long K;
    int dist_len;
    double zipf_alpha;
    char filename[100];
    // debug_level 0: only print params and final results
    //debug_level 1: print *** comparisons
    // debug_level 2: print running tally
    int debug_level;
} input_params;

int PerformRegGamma = 0;
int PerformJumping = 1;


int BloomSize;
long long Sigma;
long long K;
double zipf_alpha;

// exponential avg of jg: used to determine when converged
#define JGAVG_BOUND .001
#define JG_ALPHA 0.9





// Let phi_m(i,j) be the likelihood of going from i bits to j bits given
// msg uses m bits.  Will store this in phi_k[
double *phi_var;
double *phi_k;

// l is # of bits from cur msg so far, i, is base # bits set, j is dest # bits
// Assume 1<=l<=K, 0<=j-i<=l
int _PV(int l,int i,int j){
    if (j < i || j-i > l){
        fprintf(thefile, "ERR: PV out of scope: l: %d, i: %d, j: %d\n", l,i,j);
    }

    int idx =  i*(K+1)*(K+1)/2+TwoToOneLimited(l,j-i);
//    fprintf("%d,%d,%d maps to %d\n", l,i,j,idx);
//    fprintf("idx %d\n", idx);
    return idx;
}
// for phi_k array (i.e., where l=k)
int _PK(int i,int j){return  i*(K+1)+(j-i);}

// # for psi(l,n) array, where psi = 0 for l < n
int _PSIK(int l, int n){
    return l*(l+1)/2 + n;
}

double GenPV(int l, int i, int j){
    if (l==0){
        if (i==j){
//            fprintf("W:");
            phi_var[_PV(l,i,j)] = 1;
        }
        else {
            // fprintf("YO:\n");
            phi_var[_PV(l,i,j)] = 0;
        }
    }
    else if ((i==0 && j==0)){
        //fprintf("MO\n");
        phi_var[_PV(l,i,j)] = 0;
    }
    else {
//        fprintf("R,R,W:\n");
        phi_var[_PV(l,i,j)] = (j-i>l-1 ? 0 : phi_var[_PV(l-1,i,j)]) * (double)j /(double) BloomSize +
            (j<=i ? 0 : phi_var[_PV(l-1,i,j-1)]) * (1 - (double)(j-1) / (double)BloomSize);
/*        fprintf("pv(%d,%d,%d) = %lf * %d / B + %lf * 1-%d/B = %lf\n",
               l,i,j,phi_var[_PV(l-1,i,j)], j,
               (j<=i ? 0 : phi_var[_PV(l-1,i,j-1)]), j-1,
               phi_var[_PV(l,i,j)]);*/
        
    }
//    fprintf("R:");
    return phi_var[_PV(l,i,j)];
}

void ComputePhiK(){
    long long phi_varsize = (K+1)*(K+1)*(Sigma+K);
    long long phi_ksize = (Sigma+K)*(K+1);
    fprintf(thefile, "malloc phi_v size %lld, phi_k %lld\n", phi_varsize, phi_ksize);
    printf("malloc phi_v size %lld, phi_k %lld\n", phi_varsize, phi_ksize);
    phi_var = malloc(sizeof(double) * phi_varsize);
    phi_k = malloc(sizeof(double) * phi_ksize);
    for (int i=0; i< Sigma + K; i++){
        if (i % 1000000 == 0){
            printf("%d\n",i);
        }
        for (int l=0;l<=K;l++){
            for (int j=i;j<=i+l;j++){
                if (j < i || j-i > l){
                    fprintf(thefile, "ERR2: PV out of scope: l: %d, i: %d, j: %d\n", l,i,j);
                }
                GenPV(l,i,j);
                if (l==K){
                    phi_k[_PK(i,j)] = phi_var[_PV(K,i,j)];
                }
            }
        }
    }
    free(phi_var);
    return;
}

void FreePhik(){
    free(phi_k);
}

void VerifyP_k(){
    for (int i=0;i<Sigma+K;i++){
        double  totprob = 0;
        for (int j=i;j<=i+K;j++){
            totprob += phi_k[_PK(i,j)];
            if (phi_k[_PK(i,j)] < 0 || phi_k[_PK(i,j)]>1){
                fprintf(thefile, "Invalid phi_k(%d,%d) = %lf\n",
                       i,j, phi_k[_PK(i,j)]);
            }
        }
        if (fabs(totprob-1) > .00001){
            fprintf(thefile, "non-1 tot prob %lf for i=%d\n",totprob, i);
            for (int j=i;j<=i+K;j++){
                fprintf(thefile, "f[%d,%d] = %lf\n", i,j,phi_k[_PK(i,j)]);
            }
        }
        
    }
    return;
}
                    

double GetLucky(int bits_set){
    return phi_k[_PK(bits_set, bits_set)];
}

// steady-state probability we have i bits in the filter
double pi_i[MAXBLOOMSIZE];
void computeSteadyStateArr(){
    pi_i[0] = 1;
    double tot_pi_sum = pi_i[0]; // cumulative sum
    for (int i = 1; i <= Sigma; i++){
        double pi_guess = 0;
        for (int s = 1; s <= K; s++){
            if (i - s < 0){ // too back far in states 
                continue;
            }
            else{
                pi_guess += pi_i[i-s] * phi_k[_PK(i-s,i)];
            }
        }
        // pi_guess /= (1 - phi_k[_PK(i,i)]);
        double pi_den = (1 - phi_k[_PK(i,i)]);
        pi_i[i] = pi_guess;
        pi_i[i] = pi_i[i] / pi_den;
        tot_pi_sum += pi_i[i];
    }

    // now normalize
    double cumSum = 0;
    for (int i = 0; i <= Sigma; i++){
        pi_i[i] /= tot_pi_sum;
        // printf("pi_%d: %f\n", i, pi_i[i]);
        cumSum += pi_i[i];
    }
    // printf("%f\n",cumSum);
}

// overall false pos rate of filter for a particular M, K, Sigma
// requires computeSteadyStateArr to have run first 
double computeFalsePosRate(){
    double fpRate = 0;
    for (int i = 0; i <= Sigma; i++){
        double ep = pow((i*1.0)/BloomSize, K);
        fpRate += pi_i[i] * ep;
    }
    return fpRate;
}

// psi_k(l,n) array- probabilty n messages occupy l bits
// requires phi_k array to be computed first 
double* psi_k;
double PSI_EPSILON = 1e-6;
void ComputePsiK(){
    long long psi_ksize = ceil((Sigma+1)*(Sigma+1)/2) + Sigma + 1;
    printf("malloc psi_k size %lld\n", psi_ksize);
    psi_k = malloc(sizeof(double) * psi_ksize);
    
    // init to all zeroes 
    for (int i = 0; i < psi_ksize; i++){
        // if (i % 1000000 == 0){
        //     printf("i:%d\n",i);
        // }
        psi_k[i] = 0;
    }

    // loop through number of messages 
    for (int n = 0; n <= Sigma; n++) {
        // loop through number of bits, starting from the number of messages
        if (n % 1000 == 0){
            printf("n: %d\n",n);
        }
        double l_col_sum = 0; // for normalizing-- take sum across l for every n
        int lCutoffVal = Sigma;
        for (int l = n; l <= Sigma; l++){
            // printf("computing (%d,%d)\n", l, n);
            // base case for n = 0 
            if (n == 0){
                // psi(0,0) = 1; 0 otherwise
                if (l == 0){
                    psi_k[_PSIK(l,n)] = 1;
                }
                else{
                    psi_k[_PSIK(l,n)] = 0;
                }
                l_col_sum += psi_k[_PSIK(l,n)];
                lCutoffVal = l;
                break; // jump to next n 
            }
            // base case for l > kn
            else if (l > K * n){
                psi_k[_PSIK(l,n)] = 0;
                l_col_sum += 0;
                lCutoffVal = l;
                break; // jump to next n 
            }
            // recursive case 
            else{

                double psiSum = 0;
                
                for (int j = 1; j <= K; j++){
                    
                    if (l-j >= 0){ // we can still go back this far in bits
                        
                        psiSum += psi_k[_PSIK(l-j,n-1)] * (phi_k[_PK(l-j,l)]) / (1 - phi_k[_PK(l-j,l-j)]);
                        
                    }
                    else{ // not enough bits to go this far back 
                        break;
                    }
                }

                // printf("ind: %d\n", _PSIK(l,n));
                psi_k[_PSIK(l,n)] = psiSum;

                l_col_sum += psi_k[_PSIK(l,n)];

                
                // determine if we cut off l 
                // probability should be decreasing with l in the first place, and small
                if ((psi_k[_PSIK(l,n)] < psi_k[_PSIK(l-1,n)]) && psi_k[_PSIK(l,n)] < PSI_EPSILON){
                    lCutoffVal = l;
                    break; // jump to next n
                }
                continue;
            }
        }

        // normalize 
        // printf("n=%d\n",n);
        double normalSum = 0;
        for (int l = n; l <= lCutoffVal; l++){
            
            // very improbable cases (n close to Sigma)
            // in these cases, l_column_sum = 0 
             
            if (l_col_sum == 0) {
                if (psi_k[_PSIK(l,n)] != 0 && l != Sigma) { // something went wrong somewhere
                    printf("fatal error-- denominator in normalized phi l,n = 0!\n");
                    printf("l, n: %d, %d\n", l, n);
                    exit(0);
                }
                else{
                    // psi_k[_PSIK(Sigma,n)] = 1;
                    // normalSum = 1;
                }
            }
            else{
                psi_k[_PSIK(l,n)] /= l_col_sum;
                normalSum += psi_k[_PSIK(l,n)];
            }
            // printf("psi(%d,%d) = %f\n", l, n, psi_k[_PSIK(l,n)]);
        }
        // printf("sum: %f\n", normalSum);
        // printf("---\n");
    }
}

long long findOptimalSigma(double targetFPRate, double targetFpEpsilon, double fpr){

    long long SigmaLow = -1;
    long long SigmaHi = -1;

    while(fabs(fpr - targetFPRate) > targetFpEpsilon){
        printf("new sigma: %lld\n", Sigma);
        if (fpr > targetFPRate + targetFpEpsilon){ // too big, make smaller by reducing Sigma
            if (SigmaLow == -1){
                SigmaLow = 0;
            }
            SigmaHi = Sigma;
            Sigma = ceil((SigmaLow + SigmaHi) / 2);
        }
        else if (fpr < (targetFPRate - targetFpEpsilon)) {
            if (SigmaHi == -1) {
                SigmaHi = BloomSize;
            }
            SigmaLow = Sigma;

            Sigma = ceil((SigmaLow + SigmaHi) / 2);
        }

        if (llabs(SigmaLow - SigmaHi) <= 1){
            break;
        }

        FreePhik();
        ComputePhiK();
        VerifyP_k();
        
        computeSteadyStateArr();
        fpr = computeFalsePosRate();

    }
    return Sigma;
}

typedef struct {
    double val[2];
} SwapDouble;

int SwapActive = 0;

double Eta_i[MAX_Z];
double Eta = 1;

double JumpEta = 0;
double Eta_Jump_i[MAX_Z];

int firstEtaIter = 1;
void Init_Eta_i(){
    firstEtaIter = 1;
    Eta = 1;
    for (int i=0;i<dist_len;i++){
        Eta_i[i] = udist[i];
    }
    return;
}



// void Iter_Eta_i(){
//     Eta = 0;
//     for (int i=0;i<dist_len;i++){
// //        fprintf(thefile, "Eta[%d] before: %lf, ", i, Eta_i[i]);
//         Eta_i[i] *= (1 - udist[i]);
// //        fprintf(thefile, "after: %lf\n",  Eta_i[i]);
//         Eta += Eta_i[i];
//     }
//     return;
// }

void Iter_Eta_i(){
    if (firstEtaIter == 1){ // l = 1, Eta = 1
        firstEtaIter = 0;
        Eta = 1;
    }
    else{
        Eta = 0;
        for (int i=0;i<dist_len;i++){
            Eta_i[i] *= (1 - udist[i]);
            Eta += Eta_i[i];
        }
   }
    return;
}

void Iter_Jump_Eta_i(int M){
    JumpEta = 0;
    for (int i=0;i<dist_len;i++){
        Eta_Jump_i[i] *= pow(1 - udist[i], M);
        JumpEta += Eta_Jump_i[i];
    }
}

void deIter_Eta_i(){
    Eta = 0;
    for (int i=0;i<dist_len;i++){
        Eta_i[i] /= (1 - udist[i]);
        Eta += Eta_i[i];
    }
    return;
}

void deIter_Jump_Eta_i(int M){
    JumpEta = 0;
    for (int i=0;i<dist_len;i++){
        Eta_Jump_i[i] /= pow(1 - udist[i], M);
        JumpEta += Eta_Jump_i[i];
    }
    return;
}


double f_b[MAXBLOOMSIZE];

// For a given config, 0th entry in array will store i=sigma-K (i.e., the
// lowest possible # of entries in the frozen BF
double FrozenFilterHas_i[MAX_K];
double FrozenFilterHas_i_Jump[MAX_K];
// conditional probability that given has i bits set, jumps to or beyond sigma
// (i.e., freezes).
double FreezesFrom_i[MAX_K];

void Init_f_b(){
    f_b[0] = 1;
    for (int i=1;i<=Sigma+K;i++){
        f_b[i] = 0;
    }

    // Will also clear out Frozen
    for (int i=0;i<K;i++){
        FrozenFilterHas_i[i] = 0;
        FreezesFrom_i[i] = 1;
    // This assumes phi was already calculated
        for (int j=0;j<(K-i);j++){
            FreezesFrom_i[i] -= phi_k[_PK(Sigma-K+i,Sigma-K+i+j)];
        }
    }

    return;
}

double Alpha = 0;
double Beta = 0;
unsigned long int ell = 0;
unsigned long int lastEll = 0;


void Iter_f_b(){

    // Eta must already be bumped to ellth iteration before calling
    double new_val;
    for (int b=Sigma+K;b>=0;b--){
        new_val = (1-Eta) * f_b[b];
        for (int i=0; i<=K; i++){
            //fprintf("b-i is %d\n", b-i);
            new_val += Eta * f_b[b-i] * phi_k[_PK(b-i,b)];
            if (b-i==0){
                break;
            }
        }
        f_b[b] = new_val;
    }

    //Update frozen filter stuff
    if (1==1 || f_b[Sigma-K]>0){ // note that latter condition likely sufficient
        double frozeSum = 0;
        double overallFactor = 0;
        for (int i=0;i<K;i++){

            
            double EtaPlusOne = 0;
            for (int i=0;i<dist_len;i++){
                EtaPlusOne += (Eta_i[i] * (1 - udist[i]));
            }

            FrozenFilterHas_i[i] += f_b[i+Sigma-K] * FreezesFrom_i[i] * EtaPlusOne;
            if (FrozenFilterHas_i[i] > .00001){
                //printf("FF[%lld] at l=%ld: %.20lf\n", i+Sigma-K, ell, FrozenFilterHas_i[i]);
                //printf("f_b[%lld]: %f, Frfr: %lf\n",i+Sigma-K, f_b[i+Sigma-K], FreezesFrom_i[i]);
                frozeSum += FrozenFilterHas_i[i];
                double miss = 1 - pow(((Sigma-K+i)*1.0)/(BloomSize*1.0),K);
                overallFactor += FrozenFilterHas_i[i] * miss;
                //printf("frozeSum: %f, overallFactor: %f\n", frozeSum, overallFactor);
            }
        }
        int closeToFreeze = 0;
        for (int i = 0; i < K; i++){
            if (FrozenFilterHas_i[i] > .00001){
                 closeToFreeze = 1;
            }
        }
        if (closeToFreeze == 1){
            printf("l = %ld, eta = %f, frozeSum: %f, overallFactor: %f\n",ell, Eta, frozeSum, overallFactor);
        }

    }
/*    for (int b=0; b<= Sigma+K;b++){
        fprintf("f[%d]: %lf\n", b, f_b[b]);
        }*/
    return;
}


// Run this before recomputing f's, since depends on l-1
double PsiAndGetLucky(double *gl){ // return Psi, pass gl
    double psi = 0;
    *gl = 0;
    for (int i=0;i < Sigma; i++){
        psi += f_b[i];
        *gl += f_b[i] * GetLucky(i);
    }
    if (*gl > 1){
        fprintf(thefile, "How did Getlucky Exceed 1?\n");
    }
    return psi;
}



// The first update will be for ell = 0.
double UpdateAlphaBeta(){ // return gamma
    double gl;
    double psi;
    psi = PsiAndGetLucky(&gl);
//    fprintf(thefile, "psi: %lf, gl: %lf\n", psi, gl);
    Alpha += psi * Eta * (1 - gl);
    Beta += psi;

    // Now prep for ell++

    Iter_Eta_i();
    Iter_f_b();
    return Alpha / Beta;
}

double Xi_bound = 0;
double *first_chi;
void InitLowerValid(){
    // first_chi[i] is the probability that, if starting with Sigma-K+i bits
    // set, that the next arrival crosses the sigma threshold
    first_chi = malloc(sizeof(double)*K);
    for (int i=Sigma-K;i<Sigma;i++){
        first_chi[i-Sigma+K] =  0;
        for (int j=Sigma; j < Sigma+K; j++){
            first_chi[i-Sigma+K] += phi_k[_PK(i,j)];
        }
    }
    return;
}

double Psi_Ell(){
    double temp;
    double psi = PsiAndGetLucky(&temp);
    double psi_ell = 1 - psi;
    for (int i=Sigma-K;i<Sigma;i++){
        psi_ell += psi * f_b[i] * first_chi[i-Sigma+K];
    }
    return psi_ell;
}

int LBValid(double epsilon){
    double psi_ell = Psi_Ell();
    if (psi_ell > 1-epsilon){
        return 1;
    }
    return 0;
}

double LowerBound(){
    double temp;
    double psi = PsiAndGetLucky(&temp);
    return (Alpha+psi)/ (Beta +  psi/ Eta);
}

// Jump approximation stuff

double JumpAlpha = 0;
double JumpBeta = 0;
double jump_ell = 0;



int num_jumps = 0;


double j_b[MAXBLOOMSIZE]; // jump bit probs

// What value of M leads to an expected # of 1 new message
double last_M = 1;


void BootJumpAlphaBeta(){
    JumpAlpha = Alpha;
    JumpBeta = Beta;
    // How much "actual" ell increased 
    jump_ell = ell;
    JumpEta = Eta;
    num_jumps = 0;
    for (int i=0;i<MAXBLOOMSIZE;i++){
        j_b[i] = f_b[i];
    }
    for (int i = 0; i < dist_len; i++){
        Eta_Jump_i[i] = Eta_i[i];
    }    

    return;
}

void Iter_j_b(double short_em_is_one, int M){


    // Eta must already be bumped to ellth iteration before calling
    double new_val;
    for (int b=Sigma+K;b>=0;b--){
        new_val = short_em_is_one * j_b[b]; // instead of (1-Eta) * f_b[b];
        for (int i=0; i<=K; i++){
            //fprintf("b-i is %d\n", b-i);
            new_val += (1-short_em_is_one) * j_b[b-i] * phi_k[_PK(b-i,b)];
            if (b-i==0){
                break;
            }
        }
        j_b[b] = new_val;
    }
/*    for (int b=0; b<= Sigma+K;b++){
        fprintf("f[%d]: %lf\n", b, f_b[b]);
        }*/

    //Update frozen filter stuff
    if (1==1 || j_b[Sigma-K]>0){ // note that latter condition likely sufficient
        double frozeSum = 0;
        double overallFactor = 0;
        for (int i=0;i<K;i++){
            

            FrozenFilterHas_i[i] += j_b[i+Sigma-K] * FreezesFrom_i[i];
            if (FrozenFilterHas_i[i] > .00001){
                //printf("FF[%lld] at l=%f: %.20lf\n", i+Sigma-K, jump_ell, FrozenFilterHas_i[i]);
                //printf("f_b[%lld]: %f, Frfr: %lf\n",i+Sigma-K, f_b[i+Sigma-K], FreezesFrom_i[i]);
                frozeSum += FrozenFilterHas_i[i];
                double miss = 1 - pow(((Sigma-K+i)*1.0)/(BloomSize*1.0),K);
                overallFactor += FrozenFilterHas_i[i] * miss;
                //printf("froze at %d Sum: %f, overallFactor: %f\n", i, frozeSum, overallFactor);
            }
            // printf("M: %d\n", M);
            if (frozeSum != frozeSum){
                printf("got nan!\n");
                printf("j_b: %f, freeze from i: %f, jump eta: %f\n", j_b[i+Sigma-K], FreezesFrom_i[i], JumpEta);
                //exit(0);
            }
            
        }
        if (FrozenFilterHas_i[0] > .00001){
            printf("M: %d\n", M);
            printf("l = %f, eta = %f, frozeSum: %f, overallFactor: %f\n",jump_ell, JumpEta, frozeSum, overallFactor);
        }
    }
    return;
}

double JumpPsiAndGetLucky(double *gl){ // return Psi, pass gl
    double psi = 0;
    *gl = 0;
    for (int i=0;i < Sigma; i++){
        psi += j_b[i];
        *gl += j_b[i] * GetLucky(i);
    }
    if (*gl > 1){
        fprintf(thefile, "How did Getlucky Exceed 1?\n");
    }
    return psi;
}

double Jump_Psi_Ell(){
    if (num_jumps < Sigma / K){ // couldn't have reached yet
        return -1;
    }
    double temp;
    double psi = JumpPsiAndGetLucky(&temp);
    double psi_ell = 1 - psi;
    for (int i=Sigma-K;i<Sigma;i++){
        psi_ell += psi * j_b[i] * first_chi[i-Sigma+K];
    }
//    fprintf("Jump-psi-ell: %.20f\n", psi_ell);
    return psi_ell;
}



#define NUM_JUMP_TRACK 1000000
long unsigned int jump_track_ell[NUM_JUMP_TRACK];
double jump_track_gamma[NUM_JUMP_TRACK];



// The two following functions are another way to estimate M s.t. E[M]=1
// Currently looks like makes insignificant difference.
double New_my_eta_adjust(double new_m){
    double sum_adjusted = 0;
    for (int i=0;i<dist_len;i++){
        double my_eta_adjust = pow(1-udist[i], jump_ell+new_m);
        sum_adjusted += my_eta_adjust;

    }
    return sum_adjusted;
}

double New_variance_adjust(double new_m){
    double sum_adjusted = 0;
    for (int i=0;i<dist_len;i++){
        double my_eta_adjust = pow(pow(1-udist[i], jump_ell+new_m),2);
        sum_adjusted += my_eta_adjust;

    }
    return sum_adjusted;
}

double AdjustFunc(double e_m, double new_m){
    return e_m - New_my_eta_adjust(new_m);
//        - (e_m*e_m - New_variance_adjust(new_m));
}
        

    
        

// best_err will return how far short of E[m]=1 we are
// tries to return M such that expected # new msgs = 1
double M_Expects_1(double *best_err){
    // Start by assuming that my_eta are set to the result from the prior M
    // That is the minimum M needed this time around
    // jump_ell has already been incremented by this previous M
    double E_M = 0;
    for (int i=0;i<dist_len;i++){
        E_M += pow(1-udist[i],jump_ell);
    }
    if (E_M <= 1){
        // There is no M that can be big enough to give us E_M=1.
        // Return 0 and call it quits
        return 0;
    }
    double adjusted = 0;
    //doubling phase
    double min_M = last_M;
    double max_M = min_M/2;
    while (adjusted < 1){
        min_M = max_M;
        max_M *= 2;
        adjusted = AdjustFunc(E_M, max_M);
//        fprintf("doubling: max_M=%lf, adjusted=%lf\n", max_M, adjusted);
    }
    double best_M = min_M;
    *best_err = 1 - AdjustFunc(E_M, min_M);
    // binary searh phase
    while (fabs(*best_err) > .001){
        double cur_M = min_M + (max_M - min_M)/2;
        adjusted = AdjustFunc(E_M, cur_M);
//        fprintf("Raw EM:%lf, last-err: %lf, max: %lf, min: %lf, M: %lf, adj: %lf\n",
//E_M, *best_err, max_M, min_M, cur_M, adjusted);
        if (adjusted < 1){ // only use < 1
            if (fabs(1-adjusted)<fabs(*best_err)){
                *best_err = 1-adjusted;
                best_M = cur_M;
            }
        }
        if (adjusted > 1){
            max_M = cur_M;
        }
        else min_M = cur_M;
        if (max_M - min_M < .0001){
            fprintf(thefile, "MAX-MIN-Collision...\n");
        }
    }
    if (*best_err < 0){
        fprintf(thefile, "WARNING: best_err is negative\n");
    }
      
    return best_M;
}


double UpdateJumpAlphaBeta(){
    double gl;
    double psi;
    psi = JumpPsiAndGetLucky(&gl);
    //fprintf("E_m: %lf v 1/Eta: %le\n", E_m, 1/JumpEta);
    double short_em_is_one=0;
#ifdef _NEWVER
    // ALTERNATE DENOM:
        
    last_M = M_Expects_1(&short_em_is_one);
    JumpBeta += psi * last_M;
    
                                       
    // END ALTERNATE:
#else
    // ORIG:
    double JumpEta = 0;
    for (int i=0;i<dist_len;i++){
        JumpEta += Eta_i[i] * pow(1-udist[i], jump_ell-ell);
    }
    JumpBeta += psi * psi / JumpEta;
#endif

#ifndef _NEWVER
    // OLD:
    jump_ell += 1/JumpEta;
#else
    // New:
    jump_ell += last_M;
    //fprintf(thefile, "Adding %le to jump_ell\n", last_M);
    //fprintf("Jumped by %lf, maxM jump would be %lf\n", 1/JumpEta,last_M);
#endif
    
    // iter eta_l to eta_(l+M)
    // JumpEta = 0;
    // for (int i=0;i<dist_len;i++){
    //     JumpEta += Eta_Jump_i[i] * pow(1-udist[i], last_M);
    // } 
    JumpAlpha += psi * (1 - gl);

    // iter j(l,b) to j(l+m,b)
    Iter_j_b(short_em_is_one, last_M);
    
    // For tracking purposes now
    if (num_jumps < NUM_JUMP_TRACK){
        jump_track_ell[num_jumps] = jump_ell;
        jump_track_gamma[num_jumps] = JumpAlpha / JumpBeta;
/*        fprintf("Stored in %d: je=%d, jg=%lf\n",
               num_jumps, jump_track_ell[num_jumps],
               jump_track_gamma[num_jumps]);*/
    }
    // end tracking purposes
    num_jumps++;
    return JumpAlpha / JumpBeta;
}

double JumpLowerBound(){
    double temp;
    double psi = JumpPsiAndGetLucky(&temp);
    return (JumpAlpha + psi) / (JumpBeta + psi * last_M);
}



double HoursLeftEstimate(struct timeval start){
    if (jump_ell <= ell){
        return -1; // can't make estimate
    }
    struct timeval now;
    gettimeofday(&now,NULL);
    double elapsed = (double) (now.tv_sec - start.tv_sec);
    double est_time = (jump_ell - ell) * elapsed / ell;
    return est_time / 3600;
}
    
    




int Calculate(input_params p){
    double hrs_estimate=-1;
    double delta_estimate;
    thefile = stdout;
    BloomSize = p.BloomSize;
    Sigma = p.Sigma;
    K= p.K;
    dist_len = p.dist_len;
    zipf_alpha = p.zipf_alpha;
    int debug_level = p.debug_level;
    if (strlen(p.filename)>5){
        printf("Will use filename %s\n", p.filename);
        thefile = fopen(p.filename, "wt");
    }
    GenZipf(zipf_alpha);
    
    fprintf(thefile, "B: %d, S: %lld, K: %lld, dlen: %d, alpha: %lf\n", BloomSize, Sigma, K, dist_len, zipf_alpha);
    
    

    
    double targetFPRate = 0.0001;
    double targetFpEpsilon = targetFPRate / 10.0;

    ComputePhiK();
    VerifyP_k();

    ComputePsiK();

    printf("%f\n",psi_k[_PSIK(8,4)]);
    // exit(0);
    double sum = 0;
    for (int i = 4; i <= Sigma; i++){

        if (psi_k[_PSIK(i,4)] > .00001){
            sum += psi_k[_PSIK(i,4)];
            printf("psi_k(%d,%d)=%f\n",i,4, psi_k[_PSIK(i,4)]);
        }
    }
    printf("sum:%f\n",sum);
    exit(0);

    
    
    // compute pi_i arr-- steady-state probabilty of i bits in filter
    computeSteadyStateArr();
    // overall false positive rate for filter at default Sigma
    double fpr = computeFalsePosRate(); 
    // finds the best Sigma for the desired fpr  
    long long optSigma = findOptimalSigma(targetFPRate, targetFpEpsilon, fpr);
    fpr = computeFalsePosRate(); // at optimal sigma
    
    printf("final sigma: %lld, false pos rate: %f\n",Sigma, fpr);

    exit(0);

    InitLowerValid();
        
/*    for (int i=0;i<=K+Sigma;i++){
      for (int j=i;j<=i+K;j++){
      fprintf("phi_k(%d,%d) = %lf\n", i,j,phi_k[_PK(i,j)]);
      }
      }*/
    
    Init_Eta_i();
    Init_f_b();

    // int lMax= 1000000000;
    // int lStart = (Sigma - K) / (K); 
    // double EP_WONT_FREEZE = .00000001;
    // for(int l = 1; l <= lMax; l++){
    //     double gl;
    //     double psi = PsiAndGetLucky(&gl);
    //     if (psi < EPSILON){
    //         break;
    //     }

    //     double eta = Eta;
    //     if (eta < EP_WONT_FREEZE){
    //         printf("never freezes\n");
    //         break;
    //     }
    //     int i = Sigma - K;
    //     Iter_Eta_i();
    //     Iter_f_b();
    //     if (l < lStart){
    //         if (l % 100 == 0){
    //             printf("l:%d, psi:%f, eta:%f, lStart:%d \n", l, psi, eta, lStart);
    //         }
    //     }
    //     else{
    //         printf("l:%d, psi:%f, eta:%f \n", l, psi, eta);
    //     }
    // }
    // exit(0);
    
    struct timeval start, last, now, jump_ell_done;
    gettimeofday(&last, NULL);
    start = last;
    
    int jump_compare_idx = 0;
    int jumping=0;
    int reg_first_v = -1;
    int stop_jumping = -1; // when no longer -1, stops the jumping
    double jumpgamma;
    double gamma;
    double jump_psi_ell = 0;
    double psi_ell=0;
    int did_flip=0; // stays 0 as long as jump bound > ell bound
    double jg_avg = 1;
    while (psi_ell < 1-EPSILON){ 
        if (PerformJumping==1 && jumping==0 && Eta < ETA_START_JUMP){ 
            jumping=1; 
            BootJumpAlphaBeta(); 
        } 
        if (jumping==1 && stop_jumping < 0){ 
            jumpgamma = UpdateJumpAlphaBeta();
            jg_avg = jg_avg * JG_ALPHA + (1-JG_ALPHA) * jumpgamma;
            if (jg_avg-jumpgamma < JGAVG_BOUND * jumpgamma){
                stop_jumping = 0;
                fprintf(thefile, "STOPPING JumpGamma updates!\n");
                gettimeofday(&jump_ell_done, NULL);
            }
        } 
        if (PerformRegGamma==1 || jumping==0){ 
            gamma = UpdateAlphaBeta(); 
        } 
        gettimeofday(&now, NULL);
        unsigned long int elapsed_sec = now.tv_sec - start.tv_sec;
        psi_ell = Psi_Ell();
        char printme[1000];
        
        if (now.tv_sec > last.tv_sec || now.tv_sec - start.tv_sec < 10){
            double deltaSec = now.tv_sec - last.tv_sec;
            double deltaEll;
            if (lastEll == 0){
                deltaEll = 0;
            }
            else{
                deltaEll = ell - lastEll;
            }
            lastEll = ell;

            double lPerSec = deltaEll / deltaSec;
            last = now;
            if ((PerformRegGamma==1 || jumping==0) &&
                (debug_level >=2 || thefile != stdout)){
                
                sprintf(printme, "l: %lu, eta: %lf, gamma: %.10lf, l-ratio: %.5lf, l/sec: %f", ell, Eta, gamma, LowerBound() / gamma, lPerSec);
                printf("psi ell: %f\n", 1-psi_ell);
                if (psi_ell > 1-EPSILON){
                    sprintf(printme+strlen(printme), " (V)");
                    if (reg_first_v < 0){
                        reg_first_v = ell;
                    }
                }
                else sprintf(printme+strlen(printme), " (%lf)", 1-psi_ell);
                sprintf(printme+strlen(printme), "\n");

                if (thefile != stdout && debug_level < 2){
                    printf("%s", printme);
                }
                if (debug_level >=2){
                    fprintf(thefile, "%s",printme);
                    
                }
            }
            double prev_hr_estimate = hrs_estimate;
            hrs_estimate = HoursLeftEstimate(start);
            if (prev_hr_estimate==-1){
                delta_estimate = 0;
            }
            else delta_estimate = delta_estimate * EST_ALPHA +
                     (1 - EST_ALPHA) * (hrs_estimate - prev_hr_estimate);
            
                
            if (hrs_estimate > MAX_HOUR_EST_RUN &&
                elapsed_sec > MIN_SEC_MAKE_CALL &&
                jump_ell > MIN_JL_MAKE_CALL &&
                fabs(delta_estimate) < DELTA_MAX){
                sprintf(printme+strlen(printme),
                        "TIME-ABORT after %lu sec: Hr estimate of %.2lf too long to reach jl=%.0lf\n",
                        elapsed_sec, hrs_estimate, jump_ell);
                printf("%s",printme);
                fprintf(thefile, "%s", printme);
                return 1;
            }
            
            if (jumping==1 && (debug_level >= 2 || thefile != stdout)){
                sprintf(printme, "jl: %.1lf, nj: %d, jg: %.10f, jlr: %f", jump_ell, num_jumps, jumpgamma, jumpgamma / gamma);
                jump_psi_ell = Jump_Psi_Ell();
                if (jump_psi_ell == -1){
                    sprintf(printme+strlen(printme), " (N)");
                }
                else if (jump_psi_ell > 1-EPSILON){
                    sprintf(printme+strlen(printme), " (V)");
                }
                else sprintf(printme+strlen(printme), " (%lf)", jump_psi_ell);
                sprintf(printme+strlen(printme), "\n");

                if (hrs_estimate > 0){
                    sprintf(printme+strlen(printme),
                            "Est remaining time (hrs): %lf, delta: %lf\n", hrs_estimate, delta_estimate);
                }
                
                if (thefile != stdout && debug_level < 2){
                    printf("%s", printme);
                    printf("Cur filename is %s\n", p.filename);
                }
                if (debug_level >=2){
                    fprintf(thefile, "%s", printme);
                    
                }
                
                
            }
            
            
        }
        if (jumping==1 && PerformRegGamma==1 &&
            PerformJumping==1 && ell==jump_track_ell[jump_compare_idx]){
            if (p.debug_level >=1 || thefile != stdout){
                sprintf(printme, "*****ell:%lu, jg/g R:%lf, n-ell=%lu (%e)\n",
                       jump_track_ell[jump_compare_idx],
                       jump_track_gamma[jump_compare_idx] / gamma,
                       jump_track_ell[jump_compare_idx+1], 1-psi_ell);
            }
            if (thefile != stdout && debug_level < 1){
                printf("%s", printme);
            }
            if (debug_level >=1){
                fprintf(thefile, "%s", printme);
            }
            if (jump_track_gamma[jump_compare_idx] / gamma < 1){
                sprintf(printme+strlen(printme), "Jump-not-ub at ell=%lu\n", ell);
                did_flip = 1;
            }
            jump_compare_idx++;
        }

        ell++;
        if (jumping == 1 && (ell > jump_ell)){

            double frozeSum = 0;
            double overallFactor = 0;
            for (int i = 0; i < K; i++){
                frozeSum += FrozenFilterHas_i[i];
                double miss = 1 - pow(((Sigma-K+i)*1.0)/(BloomSize*1.0),K);
                overallFactor += FrozenFilterHas_i[i] * miss;
            }
            double gl; 
            double psi = JumpPsiAndGetLucky(&gl);
            printf("l:%ld, psi:%f, froze at i sum: %f, overall factor: %f\n",ell, psi, frozeSum, overallFactor);
            
            if (stop_jumping == 0){
                printf("jump finished and stop jumping \n");
                break;
            }
            // exit(0);
        }
    }
    fprintf(thefile, "FINAL: l: %lu, gamma: %.10lf, l-ratio: %.5lf\n",
           ell, gamma, LowerBound() / gamma);
    double jlr = jumpgamma / gamma;
    fprintf(thefile, "jl: %.1lf, nj: %d, jlr: %f\n", jump_ell, num_jumps,
           jumpgamma / gamma);

    gettimeofday(&now,NULL);
    unsigned long int elapsed = now.tv_sec - start.tv_sec;
    double e = (double) elapsed;
    if (jumping==1){
        unsigned long int jump_elapsed = jump_ell_done.tv_sec - start.tv_sec;
        double jt_ratio;
        if (jump_elapsed==0){
            jt_ratio = e * 1000000 / (double) jump_ell_done.tv_usec;
        }
        else jt_ratio = e / (double) jump_elapsed;
        fprintf(thefile, "jt_ratio: %.2f jtime: %lu sec\n", jt_ratio,
                jump_elapsed);
    }
    fprintf(thefile, "Time elapsed: %.3f hrs = %.3f min = %lu sec\n", e / 3600, e / 60, elapsed);

    if (jumping==0){
        fprintf(thefile, "GOOD: NEVER JUMPED\n");
        return 1;
    }
    if (jlr < 1){
        fprintf(thefile, "BAD by %le\n", jumpgamma/gamma-1);
        fclose(thefile);
        return 0;
    }
    else fprintf(thefile, "GOOD\n");
    if (did_flip==1){
        fprintf(thefile, "DID FLIP\n");
        fclose(thefile);
    }
    return 1;
}


