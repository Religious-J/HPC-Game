#include<iostream>
#include<stdlib.h>
#include<omp.h>
#include<cmath> 
#include"rhs.h"

#define BLOCK 64

using namespace std;

int main(int argc,char *argv[]){

 	unsigned int N1 = strtoul(argv[1], NULL, 10);
 	unsigned int N2 = strtoul(argv[2], NULL, 10);
 	unsigned int N3 = strtoul(argv[3], NULL, 10);

 	FILE *out;
    out = fopen("output.dat","w");
 	 omp_set_num_threads(8);
#pragma omp parallel num_threads(8)
{
	#pragma omp for schedule(dynamic) ordered
	for(unsigned int ii=0;ii<N1;ii+=BLOCK)
 		for(unsigned int jj=0;jj<N3;jj+=BLOCK)
 			for(unsigned int kk=0;kk<N2;kk+=BLOCK)
 				for(unsigned int i=ii;i<min(ii+BLOCK,N1);i++){
 					for(unsigned int j=jj;j<min(jj+BLOCK,N3);j++){
 						double vc=0.0;
 						for(unsigned int k=kk;k<min(kk+BLOCK,N2);k++){
 						double a,b;
 						matA(i,k,a);
 						matB(k,j,b);
 						vc += a*b;
			 			}
	 				#pragma omp ordered 
			 			fprintf(out,"%.12g\n",vc);
		 			}
	 			}
}	 
 	
    fclose(out);
    return 0;
}
