#include<iostream>
#include<stdlib.h>
#include<omp.h>
#include<cmath> 
#include"rhs.h"


using namespace std;
//
//void cau(){
//	for(int k=1;k<=10;k++){
//		for(int i=n+2;i<=(N+1)*N;i++){
//			if(i%(n+1)==0 || i%(n)==0){
//				continue;
//			}
//			x[k][i] = 0.25;
//			double re = b[i];
//			re += x[k][i-1];
//			re += x[k][i-(N+1)];
//			re += x[k-1][i+1];
//			re += x[k-1][i+N+1];
//			x[k][i]*re;
//			
//			int lop=0;
//			for(int i=n+2;i<=(N+1)*N;i++){
//				if(fabs(x[k-1][i] - x[k][i]) > e){
//					lop=1;
//					break;
//				}
//			} 
//			if(lop==0) return;
//		}
//	}
//}



int main(int argc,char *argv[]){
	
	int N = atoi(argv[1]);
 	double a[N+1][N+1];
	double b[(N+1)*(N+1)];
	double x[1000][(N+1)*(N+1)];
	
	for(int i=0;i<1000;i++)
		for(int j=0;j<(N+1)*(N+1);j++){
			x[i][j] = 0;
		}
	
	 
	double _NN = 1.0/(N*N);
	double _N = 1.0/N;
	int num=0; 
	
 	for(int i=1;i<N;i++){
 		for(int j=1;j<N;j++){
 			double xi = i*_N;
 			double yi = j*_N;
 			double value;
 			rhs(xi,yi,value);
 			b[i*(N+1)+j] = _NN * value;
		 }
	 }
	 
	  int k;
	  for(k=0;k<1000;k++){
		for(int i=N+2;i<=(N+1)*N;i++){
			if(i%(N+1)==0 || i%(N)==0){
				continue;
			}
			x[k][i] = 0.25;
			double re = b[i];
			re += x[k][i-1];
			re += x[k][i-(N+1)];
			re += x[k-1][i+1];
			re += x[k-1][i+N+1];
			x[k][i] *= re;
		}
			int lop=0;
			for(int i=N+2;i<=(N+1)*N;i++){
				if(fabs(x[k-1][i] - x[k][i]) > 1e-12){
					lop=1;
					break;
				}
			} 
			if(lop==0) break;
	  }
	
	
	FILE *out;
    out = fopen("output.dat","w");
    for(unsigned int i=0;i<(N+1)*(N+1);i++)
    	fprintf(out,"%.12g\n",x[k][i]);
    fclose(out);
    return 0;
	 
	 
	  
}
