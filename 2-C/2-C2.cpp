#include <cmath>
#include<iostream>
#include<omp.h> 
#include"rhs.h"
#include <fstream>
#include<stdlib.h>
using namespace std;


int main(int argc,char *argv[]){
	int N = atoi(argv[1]);
	double u[(N+1)*(N+1)] = {0.0};
	const int N2 = N*N;
	int flag1 = 1;
	int flag2 = 1;
	int cnt = 1;
	
	while(flag1==1 || flag2==1){
		
		flag1 = 0,flag2=0;
		#pragma omp parallel for schedule(dynamic)
		for(int len = 1; len < N;++len){
		for(int i = len; i >= 1; --i){
				int j = len - i + 1;
				double ti = 1.0*i/N;
				double tj = 1.0*j/N;
				double res;
				rhs(ti,tj,res);				
				double tmp = u[i*(N+1)+j];
				tmp=0.25*(1.0/N2*res+u[(i-1)*(N+1)+j]+u[(i+1)*(N+1)+j]+u[i*(N+1)+j-1]+u[i*(N+1)+j+1]);
				if(fabs(tmp-u[i*(N+1)+j])>1e-15){
					flag1 = 1;
				}
				u[i*(N+1)+j] = tmp;	
				
		}
	}
		
		#pragma omp parallel for schedule(dynamic)
		for(int len = N-1; len >= 2; --len){
		for(int i = N-1 ; i >= N+1-len; --i){
				int j = i - len + 2;
				double ti = 1.0*i/N;
				double tj = 1.0*j/N;
				double res;
				rhs(ti,tj,res);				
				double tmp = u[i*(N+1)+j];
				tmp=0.25*(1.0/N2*res+u[(i-1)*(N+1)+j]+u[(i+1)*(N+1)+j]+u[i*(N+1)+j-1]+u[i*(N+1)+j+1]);
				if(fabs(tmp-u[i*N+j])>1e-15){
					flag2 = 1;
				}
			u[i*(N+1)+j] = tmp;	
				
		}
		}
	
		cout << cnt++ <<endl;
		if(cnt == 100000) break;
	}
	ofstream outFile;
	outFile.open("output.dat");

	outFile.precision(12);
	outFile.setf(ios_base::showpoint);
	for(int i = 0; i < (N+1)*(N+1); ++i)
	outFile << u[i] << endl;
	outFile.close();
}
