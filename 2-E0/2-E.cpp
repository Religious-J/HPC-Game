#include<iostream>
#include<string.h>
using namespace std;

const int N = 4096*4096;
const int M = 64*64;

int main(){
//	float input[N];
//	float kernel[M];
//	float ans[N];	
	float *input = (float *)malloc(N*sizeof(float));
	float *kernel = (float *)malloc(M*sizeof(float));
	float *ans = (float *)malloc(N*sizeof(float));
	int xi,yi;
	int xk,yk;

	FILE *f1,*f2;
	f1 = fopen("input.txt","r");
	fscanf(f1,"%d",&xi);
	fscanf(f1,"%d",&yi);
	for(int i=0;i<xi*yi;i++){
		fscanf(f1,"%f",&input[i]);
	}
	f2 = fopen("weight.txt","r");
	fscanf(f2,"%d",&xk);
	fscanf(f2,"%d",&yk);
	for(int i=0;i<xk*yk;i++){
		fscanf(f2,"%f",&kernel[i]);
	}
	 fclose(f1);
	  fclose(f2);
	  
	int xa,ya;
	xa = xi-xk+1;
	ya = yi-yk+1;
	

// for(int i=0;i<xi*yi;i++){
//     cout << input[i] << endl;
// }

// for(int i=0;i<xk*yk;i++){
// 	cout << "kernel=" << kernel[i] << endl;
// }

	
	for(int i=0; i<xi-xk+1;i++){
		for(int j=0;j<yi-yk+1;j++){
			float tmp=0.0;
			for(int m=0;m<xk;m++){
				for(int n=0;n<yk;n++){
					tmp += kernel[m*yk+n] * input[(i+m)*yi+j+n];
				}
			}
			ans[i*(ya)+j] = tmp;
		}
	}
	
  
  
  
	FILE *out;
    out = fopen("output.txt","w");
    fprintf(out,"%d %d\n",xa,ya);
    for(int i=1;i<=xa*ya;i++){
    	fprintf(out,"%f ",ans[i-1]);
    	if(i % ya == 0){
    		fprintf(out,"\n");
    	}
    }
    fclose(out);
    
    free(input);
    free(kernel);
    free(ans);
} 
