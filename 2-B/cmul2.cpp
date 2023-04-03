#include<iostream>
#include<stdlib.h>
#include<omp.h>
#include<cmath> 
#include<string.h>
#include"rhs.h"

#define BLOCK 64
using namespace std;

int main(int argc,char *argv[]){

 	unsigned int N1 = strtoul(argv[1], NULL, 10);
 	unsigned int N2 = strtoul(argv[2], NULL, 10);
 	unsigned int N3 = strtoul(argv[3], NULL, 10);

	double *vc = (double *)malloc(N1*N3*sizeof(double)); 
	

	for(unsigned int i=0;i<N1*N3;i++)
		vc[i] = 0.0;
	
	//memset(vc,0.0,sizeof(double)*N1*N3);

 	
 	
 	 omp_set_num_threads(8);
#pragma omp parallel num_threads(8)
{
	#pragma omp for 

 				for(unsigned int i=0;i<N1;i++){
 					for(unsigned int k=0;k<N2;k++){
 						double a;
 						matA(i,k,a);
 						for(unsigned int j=0;j<N3;j++){
 							double b;
 							matB(k,j,b);
 							vc[i*N3+j] += a*b;
						 }
					 }
				}

}	 

 	FILE *out;
    out = fopen("output.dat","w");
    for(unsigned int i=0;i<N1*N3;i++)
    	fprintf(out,"%.12g\n",vc[i]);
    fclose(out);
    return 0;
}
