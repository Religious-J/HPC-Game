#include<iostream>
#include<stdlib.h>
#include<omp.h>
#include<cmath> 
#include<cmath> 
using namespace std;

void matA(unsigned int &i, unsigned int &j, double &value)
{
  value = std::sin(i + 2 * j);
}

void matB(unsigned int &i, unsigned int &j, double &value)
{
  value = std::cos(3 * i + j);
}


int main(int argc,char *argv[]){

 	unsigned int N1,N2,N3;
 	
 	cin >> N1 >> N2 >> N3;
 		
 	double va[N1*N2+10];
 	double vb[N2*N3+10];
 	double vc[N1*N3+10];
 	
 	for(unsigned int i=0;i<N1;i++){
 		for(unsigned int j=0;j<N2;j++){
 		double v;
 		matA(i,j,v);
 		va[i*N2+j] = v;
 		}
	 }
	 
	 for(unsigned int i=0;i<N2;i++){
 		for(unsigned int j=0;j<N3;j++){
 		double v;
 		matB(i,j,v);
 		vb[i*N3+j] = v;
 		}
	 }  

 	for(unsigned int i=0;i<N1;i++){
 		for(unsigned int j=0;j<N3;j++){
 			unsigned int index = i*N3+j;
 			vc[index] = 0;
 			for(unsigned int k=0;k<N2;k++){
 				vc[index] += va[i*N2+k] * vb[k*N3+j];
			 }
		 }
	 }
 	 for(unsigned int i=0;i<N1*N3;i++)
 	printf("%.12g\n",vc[i]);
 	
    FILE *out;
    out = fopen("output.dat","w");
    for(unsigned int i=0;i<N1*N3;i++)
    	fprintf(out,"%.12g\n",vc[i]);
    fclose(out);
    return 0;
}
