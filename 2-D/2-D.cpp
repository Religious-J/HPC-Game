#include<stdio.h>
#include<stdlib.h>
#include<string.h> 

#define data_t unsigned long long

data_t gen_A,gen_B,gen_C,n;

/*
data数组的生成的示意代码如下：
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

void data_init(){
       
}
/*
归一数的计算方式如下：
*/
data_t get_res(data_t *a,data_t n){
        for (data_t i = 0,tmp; i < n; i ++)
            for (data_t j = 0; j < n - 1 - i; j ++) {
                if (a[j] > a[j + 1]) {
                    tmp = a[j];
                    a[j] = a[j+1];
                    a[j+1] = tmp;
                }
            }
        data_t res = 0;
        for (data_t i = 0; i < n; i ++) {
            res ^= i * a[i];
        }
        return res;
}

int main(){
	
     scanf("%llu%llu%llu%llu", &gen_A, &gen_B, &gen_C, &n);
     data_t *a = (data_t *) malloc(sizeof(data_t) * n);
     for (data_t i = 0 ; i < n; i ++) {
            a[i] = gen_next();
     }
     
    data_t answer;
    answer = get_res(a,n);
    printf("%llu\n",answer);
}



#undef data_t
