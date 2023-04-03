#include<iostream>
#include<math.h>
#include<string.h>
#include<omp.h> 
#include"rhs.h"
using namespace std;

#define oover 1e-15 


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

    // initPlate(plate);
    

    //it = poissongs(plate, oover);
    double dif = 1.0f, temp;
    double _N = dif/N;
    double _NN = dif/(N*N);

    while(dif > oover){
        dif = 0.0f;             //differences at the borders are 0.0f
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
                plate[!last*(N+1)*(N+1)+i*(N+1)+j] = (_NN*res + plate[!last*(N+1)*(N+1)+(i-1)*(N+1)+j] + plate[!last*(N+1)*(N+1)+i*(N+1)+j-1] + plate[last*(N+1)*(N+1)+i*(N+1)+j+1] + plate[last*(N+1)*(N+1)+(i+1)*(N+1)+j]) * 0.25;
                temp = fabs(plate[!last*(N+1)*(N+1)+i*(N+1)+j] - plate[last*(N+1)*(N+1)+i*(N+1)+j]);
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
                plate[!last*(N+1)*(N+1)+i*(N+1)+j] = (_NN*res + plate[!last*(N+1)*(N+1)+(i-1)*(N+1)+j] + plate[last*(N+1)*(N+1)+(i+1)*(N+1)+j] + plate[!last*(N+1)*(N+1)+i*(N+1)+j-1] + plate[last*(N+1)*(N+1)+i*(N+1)+j+1] ) * 0.25;
                temp = fabs(plate[!last*(N+1)*(N+1)+i*(N+1)+j] - plate[last*(N+1)*(N+1)+i*(N+1)+j]);
                if(temp > dif) dif = temp;
            }
}
}
        last = !last; //update matrix
    }
    
    
//    for(int i=0;i<=N;i++)
//        for(int j=0;j<=N;j++){
//        cout << i << " " << j << "----" << plate[!last][i][j] << endl;
//        }
	
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
