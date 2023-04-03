#include<stdio.h> 
#include<string.h>
#include<stdlib.h>
#include<omp.h>
#include<iostream>
using namespace std; 

const int mod = 100001651;
int main(){
	int n,p;

	FILE* fi;
    fi = fopen("input.bin", "rb");
    fread(&p, sizeof(int), 1, fi);
    fread(&n, sizeof(int), 1, fi);
    int * a = (int *)malloc((n+1)*sizeof(int));
    fread(a ,sizeof(int), n, fi);
    fclose(fi);
    
    long long ans = 0;
    for(int i=0;i<n;i++){
    	ans += a[i];
    	ans = ans%mod;
	}
	ans = ans%mod;
	
	int x = ans;
    #pragma omp parallel for num_threads(p)
    for(int i=0;i<n;i++){
    	a[i]++;
	}
	printf("x = %d\n",x);
	
	FILE* fo;
    if (fo = fopen("output.bin", "wb")) {
        fwrite(&x, 1,  sizeof(int), fo);
        fwrite(a,1,sizeof(int)*n,fo); 
        fclose(fo);
    }
    free(a); 
	
}
