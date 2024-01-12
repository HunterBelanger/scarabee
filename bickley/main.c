#include <stdio.h>

extern double ki3_(double*);

int main() {

    double x = 0.;
    
    double val = ki3_(&x);

    printf("ki3(%f) = %.15f\n", x, val);

    return 0;
}