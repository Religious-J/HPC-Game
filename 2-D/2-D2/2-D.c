//omp
#include<stdio.h>
#include<stdlib.h>
#include<string.h> 
#include<math.h>
#include<omp.h> 
#include "parallel_radix_sort.h"


#define data_t unsigned long long



data_t gen_A,gen_B,gen_C,n;

/*
data閺佹壆绮嶉惃鍕晸閹存劗娈戠粈鐑樺壈娴狅絿鐖滄俊鍌欑瑓閿?
*/
data_t gen_next(){
        gen_A ^= gen_A << 31;
        gen_A ^= gen_A >> 17;
        gen_B ^= gen_B << 13; 
        gen_B ^= gen_B >> 5;
        gen_C ++;
        gen_A ^=gen_B;
        gen_B ^= gen_C;
        return gen_A;        
}

/*
瑜版帊绔撮弫鎵畱鐠侊紕鐣婚弬鐟扮础婵″倷绗呴敍?
*/


int main(){
    
     scanf("%llu%llu%llu%llu", &gen_A, &gen_B, &gen_C, &n);
     data_t *a = (data_t *) malloc(sizeof(data_t) * n);

     for (data_t i = 0 ; i < n; i ++) {
            a[i] = gen_next();
            //printf("%llu\n",a[i]);
     }
     
     
     // sort(a,0,n);
     parallel_radix_sort::SortKeys(a, n);
     
        data_t res = 0;
#pragma omp parallel for reduction(^:res)
        for (data_t i = 0; i < n; i ++) {
            res ^= i * a[i];
        }

    printf("%llu\n",res);
    
    free(a); 
}



#undef data_t
