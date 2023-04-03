## PART 1 - 高性能计算简介

### B. 实验室的新机器

```shell
#!/bin/bash
srun ./program $1 > output.dat
seff $(cat job_id.dat) > seff.dat
```

### C. 小北问答：超速版

```yaml
Q1: Henri
Q2: 1987
Q3:
  B1: 多
  B2: 强
Q4:
  B1: 低
  B2: 弱
Q5:
  B1: N*N*N
  B2: N*N*N
Q6:
  B1: 3*N*N
  B2: N*N
Q7:
  B1: N*N*N+2*N*N
  B2: N*N
Q8: 神威太湖之光
Q9:
  B1: A
  B2: B
  B3: D
Q10:
  B1: 2016
  B2: 杨超
  B3: 北京大学
```

### D. 简单题

#### code 1

使用 mmap

```c++
#include <iostream>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <cstring>
#include <thread>
#include <sys/stat.h>
using namespace std;
constexpr int mod = 100001651;
int main() {
    unsigned int P, N;
    int fd = open("input.bin", O_RDWR);
    struct stat sb;
    fstat(fd, &sb);
    unsigned int *data = (unsigned int*) mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    P = data[0];
    N = data[1];
    long long sum = 0;
    for (unsigned int i = 0; i < N; i++) {
        sum += data[i + 2];
        sum %= mod;
    }
    #pragma omp parallel for
    for (int i = 0;i<N; i++) {
        data[i+2]++;
    }
    int fd2 = open("output.bin", O_RDWR | O_CREAT | O_TRUNC, 0666);
    ftruncate(fd2, (N + 1) * sizeof(int));
    int *output = (int*) mmap(NULL, (N + 1) * sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, fd2, 0);
    output[0] = sum;
    memcpy(output + 1, data + 2, N * sizeof(unsigned int));
    munmap(output, (N + 1) * sizeof(unsigned int));
    munmap(data, (N + 3) * sizeof(unsigned int));
    close(fd);
	close(fd2);
    return 0;
}  
```

#### code 2

使用 fread/fwrite

```c++
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
```

#### code 3

使用 ifstream/ofstream

```c++
#include <stdio.h>
#include <stdint.h>     // included for uint8_t
#include <iostream>
#include <fstream>
#include<vector>
using namespace std;
int main(int argc, char const *argv[])
{		
	vector<int> c;
	ifstream ifs("input.bin",ios_base::in|ios_base::binary);
	int geoCount[2];
	ifs.read((char*)&geoCount, sizeof(int)*2);
	int p = geoCount[0];
	int n = geoCount[1];
	printf("%d,%d\n",p,n);
    
 	long long sum = 0;
 	for(int i = 0; i < n;++i){
 		int binary_tmp;
 		ifs.read((char*)&binary_tmp, sizeof(int));
 		sum += binary_tmp;
 		if(sum > 100001651) sum%=100001651;
 		c.push_back(binary_tmp+1);
	 }
	sum %= 100001651;
	int res = int(sum);
    
	ofstream out("output.bin",ios::binary);
	out.write((char*)&res,sizeof(int));
	for(int i = 0; i < n; ++i){
		int binary_tmp = c[i];
		out.write((char*)&binary_tmp,sizeof(int));
	}
	out.close();
    return 0;
}
```



## PART 2 - 并行与大规模

### A. 求积分！

```c++
#include<iostream>
#include<stdlib.h>
#include<omp.h> 
#include"rhs.h"
using namespace std;
int main(int argc,char *argv[]){
    int N = atoi(argv[1]);
    double dif = 1.0; 
    double NN = dif/N;
    double answer=0.0;
    omp_set_num_threads(8);
#pragma omp parallel num_threads(8)
{
	#pragma omp for reduction(+:answer)
    for(int i=1;i<=N;i++){
        double xi = 1.0*(i-0.5)*NN;
        double value; 
        rhs(xi,value);
        answer += value;
    }
} 
    answer = answer*NN;
    FILE *out;
    out = fopen("output.dat","w");
    fprintf(out,"%.12g\n",answer);
    fclose(out);
    return 0;
}
```

### B. 乘一乘！

```c++
#include <iostream>
#include <fstream>
#include <omp.h>
#include "rhs.h"
using namespace std;
int main(int argc, char* argv[]) {
    unsigned int N1 = strtoul(argv[1], NULL, 10);
 	unsigned int N2 = strtoul(argv[2], NULL, 10);
 	unsigned int N3 = strtoul(argv[3], NULL, 10);
    double *A = (double*)malloc(sizeof(double) * (N1 * N2)+10);
    double *B = (double*)malloc(sizeof(double) * (N3 * N2)+10);
    double *C = (double*)malloc(sizeof(double) * (N1 * N3)+10);
    for (unsigned int i=0;i<N1;i++) {
        for (unsigned int t=0;t<N2;t++) {
            matA(i, t, A[i * N2 + t]);
        }
    }
    for (unsigned int i=0;i<N2;i++) {
        for (unsigned int t=0;t<N3;t++) {
            matB(i, t, B[i * N3 + t]);
        }
    }
    
    #pragma omp parallel for
    for (unsigned int i = 0; i < N1; i++) {
        for (unsigned int j = 0; j < N3; j++) {
            for (unsigned int k = 0; k < N2; k++) {
            	C[i * N3 + j]+=A[i * N2 + k] * B[k * N3 + j];
            }
        }
    }
    
    FILE *fp = fopen("output.dat", "w");
    for (unsigned int i = 0; i < N1; i++) {
        for (unsigned int j = 0; j < N3; j++) {
			fprintf(fp, "%.11e\n", C[i * N3 + j]);
        }
    }
    fclose(fp);
    return 0;
}

```

### C. 解方程！

#### code 1

```c++
#include<iostream>
#include<math.h>
#include<string.h>
#include<omp.h> 
#include"rhs.h"
using namespace std;
#define epic 1e-15 
int main(int argc, char *argv[]){
    int N = atoi(argv[1]);
    double *plate = (double *)malloc(2*(N+1)*(N+1)*sizeof(double));
    omp_set_num_threads(8);
#pragma omp parallel num_threads(8)
{
	#pragma omp for
    for(int i=0;i<=2*(N+1)*(N+1);i++){
		plate[i] = 0.0;
	}
}	 
    int last=0;
    //initPlate(plate);
    //it = poissongs(plate, oover);
    double dif = 1.0, temp;
    double _N = dif/N;
    double _NN = dif/(N*N);
    while(dif > epic){
        dif = 0.0;             //differences at the borders are 0.0f
#pragma omp parallel sections
{
#pragma omp section
{
		//Black
        for(int i = 1; i < N; i ++)
            for(int j = 2-i%2; j < N; j +=2){
                double res;
                double xi = 1.0*i*_N;
                double yj = 1.0*j*_N;
                rhs(xi,yj,res);
                plate[!last*(N+1)*(N+1)+i*(N+1)+j] = (_NN*res + plate[last*(N+1)*(N+1)+(i-1)*(N+1)+j] + plate[last*(N+1)*(N+1)+(i+1)*(N+1)+j] + plate[last*(N+1)*(N+1)+i*(N+1)+j-1] + plate[last*(N+1)*(N+1)+i*(N+1)+j+1] ) * 0.25;
                temp = plate[!last*(N+1)*(N+1)+i*(N+1)+j] - plate[last*(N+1)*(N+1)+i*(N+1)+j];
                if(temp < 0) temp *= -1; 
                if(temp > dif) dif = temp;
            }
}
#pragma omp section
{
		//red
        for(int i = 1; i < N; i ++)
            for(int j = 1+i%2; j < N; j +=2){
                double res;
                double xi = 1.0*i*_N;
                double yj = 1.0*j*_N;
                rhs(xi,yj,res);
                plate[!last*(N+1)*(N+1)+i*(N+1)+j] = (_NN*res + plate[!last*(N+1)*(N+1)+(i-1)*(N+1)+j] + plate[!last*(N+1)*(N+1)+(i+1)*(N+1)+j] + plate[!last*(N+1)*(N+1)+i*(N+1)+j-1] + plate[!last*(N+1)*(N+1)+i*(N+1)+j+1] ) * 0.25;
                temp = plate[!last*(N+1)*(N+1)+i*(N+1)+j] - plate[last*(N+1)*(N+1)+i*(N+1)+j];
                if(temp < 0) temp *= -1; 
                if(temp > dif) dif = temp;
            }
}
}
        last = !last; //update matrix
    }
    FILE *out;
    out = fopen("output.dat","w");
    for(int i=0;i<=N;i++)
        for(int j=0;j<=N;j++){
            fprintf(out,"%.12g\n",plate[!last*(N+1)*(N+1)+i*(N+1)+j]);
        }
    fclose(out);
    free(plate);
    return 0;
}
```

#### code 2

```c++
#include <iostream>
#include <fstream>
#include <omp.h>
#include <cmath>
#include "rhs.h"
using namespace std;
const double EPS=1e-15;
int main(int argc, char* argv[]) {
    int N;
    N = atoi(argv[1]);
    double U[N+5][N+5],f[N+5][N+5];
    for(int i=0;i<=N;i++)
    	for(int j=0;j<=N;j++){
    		U[i][j]=0;
		}
	for(int i=1;i<N;i++)
    	for(int j=1;j<N;j++){
    		double x=1.0*i/N,y=1.0*j/N;
    		double va;
    		rhs(x,y,va);
    		f[i][j]=va/N/N;
		}
    while(1){
    	double eps=0;
    	#pragma omp parallel for
    	for(int i=1;i<N;i++){
    		for(int j=1;j<N;j++){
    			double temp=U[i][j];
    			U[i][j]=(U[i-1][j]+U[i+1][j]+U[i][j+1]+U[i][j-1]+f[i][j])*0.25;
    			if (eps<fabs(U[i][j]-temp)) {
                    #pragma omp critical
                    eps=fabs(U[i][j]-temp);
                }
			}
		}	
		if(eps<EPS) break;
	}
    FILE *fp = fopen("output.dat", "w");
    for (int i = 0;i<=N;i++) {
        for (int j = 0;j<=N;j++) {
			fprintf(fp, "%.12lf\n", U[i][j]);
        }
    }
    fclose(fp);
    return 0;
}

```

### D. 道生一

```shell
#!/bin/bash
export OMP_PROC_BIND=TRUE
g++ -std=c++11 -Ofast  -fopenmp -mfma  -ftree-vectorize -ffast-math  -funroll-all-loops -mavx2 -mtune=native -march=native 2-D.c -o answer
```

```c++
//omp
#include<stdio.h>
#include<stdlib.h>
#include<string.h> 
#include<math.h>
#include "parallel_radix_sort.h"

#define data_t unsigned long long
data_t gen_A,gen_B,gen_C,n;
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
```

```c++
// parallel_radix_sort.h
#ifndef PARALLEL_RADIX_SORT_H_
#define PARALLEL_RADIX_SORT_H_
......
https://github.com/iwiwi/parallel-radix-sort/blob/master/parallel_radix_sort.h
```

### E. 卷？寄！

```c++
#include<iostream>
#include<string.h>
#include<omp.h>
#include <immintrin.h>  // AVX512
using namespace std;
const int N = 4096*4096;
const int M = 64*64;
int main(){
	float *input = (float *)aligned_alloc(64, N*sizeof(float));
	float *kernel = (float *)aligned_alloc(64, M*sizeof(float));
    float *ans = (float *)aligned_alloc(64, N*sizeof(float));
    //  printf("default-aligned addr:   %p\n", (void*)input);
    //  printf("default-aligned addr:   %p\n", &input[1]);
    //  printf("default-aligned addr:   %p\n", (void*)kernel);
    //  printf("default-aligned addr:   %p\n", (void*)ans);
	int xi,yi,xk,yk;
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

	for(int i=0; i<xi-xk+1; i++){
		for(int j=0; j<yi-yk+1; j++){
         __m512 tmp = _mm512_setzero_ps();
			for(int m=0; m<xk; m++){
				for(int n=0; n<yk; n+=16){
                    tmp = _mm512_add_ps(tmp,_mm512_mul_ps(_mm512_loadu_ps(&kernel[m*yk+n]),_mm512_loadu_ps(&input[(i+m)*yi+j+n])));
                    // or
           			// tmp = _mm512_fmadd_ps(_mm512_loadu_ps(&kernel[m*yk+n]),_mm512_loadu_ps(&input[(i+m)*yi+j+n]));
				}
			}
            ans[i*(ya)+j] = _mm512_reduce_add_ps(tmp);
        	// or
            // ans[i*(ya)+j] = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+tmp[6]+tmp[7]+tmp[8]+tmp[9]+tmp[10]+tmp[11]+tmp[12]+tmp[13]+tmp[14]+tmp[15];
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
```

### F. MPI算个PI

#### code 1

```c++
#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<math.h>
double f(double x) {                                    // inline
	double y = 2 * x + 1;
	double z = pow(-1, x);
	double h1 = 4.0;
	double h2 = 5.0;
	double h3 = 239.0;
	return h1 * z / y*(h1 / pow(h2, y) - 1 / pow(h3, y));
}
int main(int argc, char* argv[])
{
	int myid, numprocs, namelen;
	double pi, sum, x, *temp;
	int n = 20000;                                      // n only 10000
	char processor_name[MPI_MAX_PROCESSOR_NAME];
    
	MPI_Init(&argc, &argv);       					 	// starts MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);  			 	// get current process id
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);        	// get number of processes
	MPI_Get_processor_name(processor_name, &namelen);

	if (myid == 0) {
		temp = (double*)malloc(sizeof(double)*numprocs);
	}
    
	MPI_Bcast(&n, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD); 
	sum = 0.0, pi = 0.0;
	for (long long i = myid; i <= n; i += numprocs) {
		sum += f(i);
	}
	MPI_Gather(&sum, sizeof(sum), MPI_BYTE, temp, sizeof(sum), MPI_BYTE, 0, MPI_COMM_WORLD);
    
	if (myid == 0) {
		for (int i = 0; i < numprocs; i++) {
			pi += temp[i];
		}
		FILE *out;
    	out = fopen("output.txt","w");
    	fprintf(out,"%.15g\n",pi);
    	fclose(out);
		free(temp);
	}
	MPI_Finalize();
	return 0;
}
```

#### code 2

```c++
#include<stdio.h>
#include<mpi.h>
int main(int argc,char *argv[])
{
    int myid,np;
    long long int i,j;
    double pi=0.0; 
    double fVal;								//fVal代表取Xi所对应的函数值   4/(1+x^2) 即每个矩形的高度
    long long int n = 6000000000;
    MPI_Status status;
    double fl = 1.0;
    double h=fl/n; 								//每个矩形的宽度
    double local=0.0;							//每个进程计算的面积和 
    double start,end,xi;
    MPI_Init(&argc,&argv);   
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    for(i=myid;i<n;i+=np) 	  					//利用np个进程同时计算各部分矩形面积
    {
            xi=(i+0.5)*h;
            fVal=4.0/(1.0+xi*xi);				//得到f(xi) 
            local+=fVal ; 
    }
    local=local*h;								//得到该进程所计算的面积 
    //进程号！=0的进程计算的结果发送到进程0上面 
    if(myid!=0)
    {   
        MPI_Send(&local,1,MPI_DOUBLE,0,myid,MPI_COMM_WORLD); 
    }
    if(myid==0) 								//进程号为0就累加每个进程的计算结果 
    {
        pi=local;								//得到进程0的值 后面接收就会覆盖这个值 
        for(j=1;j<np;j++)
           {
                MPI_Recv(&local,1,MPI_DOUBLE,j,j,MPI_COMM_WORLD,&status); //把其他进程的结果发送到local中 
                pi+=local;						//得到所有的面积和 
           }
    }
    if(myid==0)
    {
        FILE *out;
    	out = fopen("output.txt","w");
    	fprintf(out,"%.14f\n",pi);
    	fclose(out);
    }
    MPI_Finalize();
    return 0;
} 
```

