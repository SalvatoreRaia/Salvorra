#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    int i,N,seed;
    printf("FEDERICO SACCO");
    scanf("%d", &N);
    float T[N], h[N], r1[N], r2[N], s1 = 0.05, Tt[N], p[N], pt = 1000, ht[N];
    Tt[0]=300;
    ht[0]=20;
   for(i=1;i<N;i++){
        Tt[i] = (300 + i*20);
        ht[i] = (Tt[i]/(Tt[i-1]))*(ht[i-1]);
    }
    for(i=0; i<N; i++){
        r1[i] = (float)rand()/RAND_MAX;
        r2[i] = (float)rand()/RAND_MAX;
        T[i] = sqrt(-2*log(r1[i]))*cos(2*M_PI*r2[i]);
        h[i] = sqrt(-2*log(r1[i]))*sin(2*M_PI*r2[i]);
        T[i]=Tt[i]+s1*Tt[i]*T[i];
        h[i]=ht[i]+h[i];
    }
    printf("SUCA");
}

