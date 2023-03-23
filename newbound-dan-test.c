#include "newbound-dan.c"

int main(int argc, char **argv){ // params: M, sigma, k, zipf N, zipf s
    input_params p;
    printf("argc: %d\n", argc);

    const int mArrLen = 1;
    const int sArrLen = 1;
    int mArr[mArrLen] = {500};
    float sArr[sArrLen] = {1};
    for (int j = 0; j < sArrLen; j++){
        p.BloomSize = mArr[0];
        p.Sigma = atoll(argv[2]);
        p.K = atoll(argv[3]);
        p.dist_len = atoi(argv[4]);
        p.zipf_alpha = sArr[j];
        if (argc>6){
            p.debug_level = atoi(argv[6]);
        }
        else p.debug_level = 2;
        if (argc>7){
            strcpy(p.filename, argv[7]);
        }
        Calculate(p);
    }
    
    return 1;   
}