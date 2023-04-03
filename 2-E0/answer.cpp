

#include<iostream>
#include<string.h>
#include<omp.h>
#include <immintrin.h>  // avx, 16¸ö256Î»¼Ä´æÆ÷

void* aligned_malloc(size_t required_bytes, size_t alignment)
{
    int offset = alignment - 1 + sizeof(void*);
    void* p1 = (void*)malloc(required_bytes + offset);
    if (p1 == NULL)
        return NULL;
    void** p2 = (void**)( ( (size_t)p1 + offset ) & ~(alignment - 1) );
    p2[-1] = p1;
    return p2;
}

void aligned_free(void *p2)
{
    void* p1 = ((void**)p2)[-1];
    free(p1);
}


using namespace std;

const int N = 4096*4096;
const int M = 64*64;

int main(){
	
//	float input[N];
//	float kernel[M];
//	float ans[N];
	float *input = (float *)aligned_alloc(64, N*sizeof(float));
	float *kernel = (float *)aligned_alloc(64, M*sizeof(float));
    	float *ans = (float *)aligned_alloc(64, N*sizeof(float));
    	
    	
    	//  printf("default-aligned addr:   %p\n", (void*)input);
    	//  printf("default-aligned addr:   %p\n", &input[1]);
    	//  printf("default-aligned addr:   %p\n", (void*)kernel);
    	//  printf("default-aligned addr:   %p\n", (void*)ans);
    	
	// float *input = (float *)_mm_malloc(N*sizeof(float),64);
	// float *kernel = (float *)_mm_malloc(M*sizeof(float),64);
	// float *ans = (float *)_mm_malloc(N*sizeof(float),64);
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
		//float temp = 0.0;
        //__m512 tmp = _mm512_loadu_ps(ans+i*(ya)+j);
         __m512 tmp = _mm512_setzero_ps();
			for(int m=0;m<xk;m++){
				for(int n=0;n<yk;n+=16){
					//tmp = _mm512_dp_ps(_mm512_load_ps(&kernel[m*yk+n]),_mm512_load_ps(&input[(i+m)*yi+j+n]),0xf1);
                    tmp = _mm512_add_ps(tmp,_mm512_mul_ps(_mm512_loadu_ps(&kernel[m*yk+n]),_mm512_loadu_ps(&input[(i+m)*yi+j+n])));
					//temp += kernel[m*yk+n] * input[(i+m)*yi+j+n]; 
					 //temp += kernel[m*yk+n] * input[(i+m)*yi+j+n];
					// __m512 tmp = _mm512_loadu_ps(&kernel[m*yk+n]);
					 //cout << kernel[m*yk+n] << " " << kernel[m*yk+n+1] << " " << kernel[m*yk+n+2] << endl;
					 //cout << "loadu: " << tmp[0] <<" "<< tmp[1] <<" " <<tmp[2] << endl;
				}
			}
	//		cout << "tmp: " << tmp[0] << " " << tmp[1] << " " << tmp[2] <<  " " << tmp[3] << endl;
        //_mm512_storeu_ps(ans+i*(ya)+j,tmp);
       
        	ans[i*(ya)+j] = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+tmp[6]+tmp[7]+tmp[8]+tmp[9]+tmp[10]+tmp[11]+tmp[12]+tmp[13]+tmp[14]+tmp[15];
        //cout << ans[i*(ya)+j] << endl;
			//ans[i*(ya)+j] = temp;
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
