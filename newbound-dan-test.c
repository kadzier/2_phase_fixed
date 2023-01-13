#include "newbound-dan.c"

int main(int argc, char **argv){
    input_params p;
    printf("argc: %d\n", argc);
    p.BloomSize = atoi(argv[1]);
    p.Sigma = atoll(argv[2]);
    p.K = atoll(argv[3]);
    p.dist_len = atoi(argv[4]);
    p.zipf_alpha = atof(argv[5]);
    if (argc>6){
        p.debug_level = atoi(argv[6]);
    }
    else p.debug_level = 2;
    if (argc>7){
        strcpy(p.filename, argv[7]);
    }
    Calculate(p);
    return 1;
    
}
