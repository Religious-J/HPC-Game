#include<iostream>
#include<stdlib.h>
#include<omp.h> 
#include"rhs.h"
using namespace std;

int main(int argc,char *argv[]){
    int N = atoi(argv[1]);
    double dif = 1.0; 
    double NN = dif/N;
    double answer=0.0;
    omp_set_num_threads(8);
    #pragma omp parallel num_threads(8)
	{
	#pragma omp for reduction(+:answer)
    for(int i=1;i<=N;i++){
        double xi = 1.0*(i-0.5)*NN;
        double value; 
        rhs(xi,value);
        answer += value;
    }
    } 
    answer = answer*NN;
    FILE *out;
    out = fopen("output.dat","w");
    fprintf(out,"%.12g\n",answer);
    fclose(out);
    return 0;
}
