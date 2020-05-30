#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float *NR (float xn, float e) {
    float xn1, a= 0;
    int k = 0;
    do {
        xn1 = xn - (pow(xn,2)+xn-5)/(2*xn +1);
        a = fabs(xn1 - xn);
        k+=1 ;
        printf("x_n:%f\n\n Iterazione:%d\n\n", xn1, k);
        xn = xn1;
    } while (a>e);
    printf("eps :%f\n\n", a);
    static float result[2];
    result[0] = xn1;
    result[1] = a;
    return result;
}

int main() {
    float xn,eps, fxn = 0;
    float *pnt;
    pnt = NR(1, 0.010);
    xn = *(pnt + 0);
    eps = *(pnt + 1);
    printf("xn : %f +- %f", xn, eps);

}
