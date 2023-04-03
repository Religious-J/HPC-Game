#include<stdio.h>
#include<stdlib.h>
#include<string.h> 
#include<mpi.h> 
#define data_t unsigned long long

data_t gen_A,gen_B,gen_C,n;


/********** Merge Function **********/
void merge(data_t *a, data_t*b, data_t l, data_t m, data_t r) {
	
	int h, i, j, k;
	h = l;
	i = l;
	j = m + 1;
	
	while((h <= m) && (j <= r)) {
		
		if(a[h] <= a[j]) {
			
			b[i] = a[h];
			h++;
			
			}
			
		else {
			
			b[i] = a[j];
			j++;
			
			}
			
		i++;
		
		}
		
	if(m < h) {
		
		for(k = j; k <= r; k++) {
			
			b[i] = a[k];
			i++;
			
			}
			
		}
		
	else {
		
		for(k = h; k <= m; k++) {
			
			b[i] = a[k];
			i++;
			
			}
			
		}
		
	for(k = l; k <= r; k++) {
		
		a[k] = b[k];
		
		}
		
	}

/********** Recursive Merge Function **********/
void mergeSort(data_t *a, data_t *b, data_t l, data_t r) {
	
	int m;
	
	if(l < r) {
		
		m = (l + r)/2;
		
		mergeSort(a, b, l, m);
		mergeSort(a, b, (m + 1), r);
		merge(a, b, l, m, r);
		
		}
		
	}

/*
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



int main(int argc, char *argv[]){
	
     scanf("%llu%llu%llu%llu", &gen_A, &gen_B, &gen_C, &n);
     data_t *a = (data_t *) malloc(sizeof(data_t) * n);
     for (data_t i = 0 ; i < n; i ++) {
            a[i] = gen_next();
     }
     
 	

	/********** Initialize MPI **********/
	int world_rank;
	int world_size;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		
	/********** Divide the array in equal-sized chunks **********/
	data_t  size = n/world_size;
	
	/********** Send each subarray to each process **********/
	data_t  *sub_array = malloc(size * sizeof(data_t ));
	MPI_Scatter(a, size, MPI_UNSIGNED_LONG, sub_array, size, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	
	/********** Perform the mergesort on each process **********/
	data_t  *tmp_array = malloc(size * sizeof(data_t ));
	mergeSort(sub_array, tmp_array, 0, (size - 1));
	
	/********** Gather the sorted subarrays into one **********/
	data_t  *sorted = NULL;
	if(world_rank == 0) {
		
		sorted = malloc(n * sizeof(data_t));
		
		}
	
	MPI_Gather(sub_array, size,MPI_UNSIGNED_LONG, sorted, size, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	
	/********** Make the final mergeSort call **********/
	if(world_rank == 0) {
		
		data_t  *other_array = malloc(n * sizeof(data_t));
		mergeSort(sorted, other_array, 0, (n - 1));
		
			
		/********** Clean up root **********/
		
		free(other_array);
			
		}
	
	/********** Clean up rest **********/
	free(a);
	free(sub_array);
	free(tmp_array);
	
	/********** Finalize MPI **********/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    
        data_t res = 0;
        for (data_t i = 0; i < n; i ++) {
            res ^= i * sorted[i];
        }
       
    printf("%llu\n",res);
}



#undef data_t
