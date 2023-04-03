//omp
#include<stdio.h>
#include<stdlib.h>
#include<string.h> 
#include<math.h>
#include<omp.h>
#define data_t unsigned long long

data_t Partition(data_t* data, int start, int end)   //鍒掑垎鏁版嵁
{
    data_t temp = data[start];   //浠ョ涓€涓厓绱犱负鍩哄噯
    while (start < end) {
        while (start < end && data[end] >= temp)end--;   //鎵惧埌绗竴涓瘮鍩哄噯灏忕殑鏁?
        data[start] = data[end];
        while (start < end && data[start] <= temp)start++;    //鎵惧埌绗竴涓瘮鍩哄噯澶х殑鏁?
        data[end] = data[start];
    }
    data[start] = temp;   //浠ュ熀鍑嗕綔涓哄垎鐣岀嚎
    return start;
}

void quickSort(data_t* data, int start, int end)  //并行快排
{
    if (start < end) {
        data_t pos = Partition(data, start, end);
        #pragma omp parallel sections    //设置并行区域
        {
            #pragma omp section          //该区域对前部分数据进行排序
            quickSort(data, start, pos - 1);
            #pragma omp section          //该区域对后部分数据进行排序
            quickSort(data, pos + 1, end);
        }
    }
}


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
data_t get_res(data_t *a,data_t n){

//         for (data_t i = 0; i < n; i ++){
//             data_t tmp;
//             for (data_t j = 0; j < n - 1 - i; j ++) {
//                 if (a[j] > a[j + 1]) {
//                     tmp = a[j];
//                     a[j] = a[j+1];
//                     a[j+1] = tmp;
//                 }
//             }
//         }

        //quick_sort(a,0,n-1);
        
        omp_set_num_threads(2);   //璁剧疆绾跨▼鏁?
    	quickSort(a , 0, n-1);
        

        data_t res = 0;
        int num = omp_get_max_threads();
        printf("%d\n",num);
	omp_set_num_threads(num); 
#pragma omp parallel for reduction(^:res)
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
            //printf("%llu\n",a[i]);
     }
     
    data_t answer;
    answer = get_res(a,n);
    printf("%llu\n",answer);
}



#undef data_t
