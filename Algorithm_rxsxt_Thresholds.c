#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <time.h>

int additions_tensor_rxs(int r, int s, int K, double *P[K][3]);
int additions_tensor_sxt(int s, int t, int K, double *P[K][3]);
int additions_tensor_rxt(int r, int t, int K, double *P[K][3]);
double f_ord(int r, int s, int t, double l);
double f_op(int r, int s, int t, double l, int K, int a, int b, int c);
int min_integer(int r, int s, int t, int K, int a, int b, int c, double *P[K][3]);


void matrix_of_zeros(int m, int n, double **M);
void random_matrix(int m, int n,double **M);
void allocate_matrix_memory(int m, int n, double ***M);
void print_matrix(int m, int n, double **M);
void ordinary_product(int m, int n, int p, double **A, double **B, double **C);
void matrices_combination(int m, int n, double **A, int f, int c, double **A1, double *coef);
void matrix_recombination(int m, int n, double **C, int f, int c, double **A1, double *coef);
void algorithm_rxsxt_normal(int m, int n, int p, double **A, double **B, double **C,int r, int s, int t, int K, double *P[K][3]);
void algorithm_rxsxt_optim(int m, int n, int p, double **A, double **B, double **C,int r, int s, int t, int K, double *P[K][3], int L0);
void free_matrix(int m, double **M);


int min_int;
int main(){
    srand (time(NULL));
    ///Llegir tensor
     char tensor_file[20];

    printf("Dimensions r s t: ");
    int r, s, t, K;
    scanf("%d %d %d", &r, &s, &t);
    sprintf(tensor_file, "tensor_%dx%dx%d.txt", r, s, t);

    FILE *tensor;
    tensor = fopen(tensor_file, "r");

    if(tensor == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    fscanf(tensor,"%d",&K);

    double *P[K][3];

    for(int i=0;i<K;i++){
        P[i][0]=(double*)malloc((r*s)*sizeof(double));
        P[i][1]=(double*)malloc((s*t)*sizeof(double));
        P[i][2]=(double*)malloc((r*t)*sizeof(double));
        if(P[i][0] == NULL || P[i][1] == NULL || P[i][2] == NULL){
            printf("Memory allocation failed\n");
            return 1;
        }
    }

    for (int i=0;i<K;i++) {
        for (int k=0;k<r*s;k++)
            fscanf(tensor, "%lf", &P[i][0][k]);
        for(int k=0;k<s*t;k++)
            fscanf(tensor, "%lf", &P[i][1][k]);
        for (int k=0;k<r*t;k++)
            fscanf(tensor, "%lf", &P[i][2][k]);
    }

    fclose(tensor);

    clock_t start, stop;
    srand (time(NULL));
    printf("TENSOR %dx%dx%d\n",r,s,t);

    int a,b,c,L0;
    a=additions_tensor_rxs(r,s,K,P);
    b=additions_tensor_sxt(s,t,K,P);
    c=additions_tensor_rxt(r,t,K,P);
    L0=min_integer(r,s,t,K,a,b,c,P); //Global variable
    printf("L0: %d\n",L0);
    int l;
    for(l=2;l<10;l++){
    printf("l= %d\n",l);
        int m=pow(r,l);
        int n=pow(s,l);
        int p=pow(t,l);

        for(int i=0;i<4;i++){
            double **A;
            allocate_matrix_memory(m,n,&A);
            random_matrix(m,n,A);
            //print_matrix(m,n,A);
            //printf("\n\n");

            double **B;
            allocate_matrix_memory(n,p,&B);
            random_matrix(n,p,B);
            //print_matrix(m,n,B);
            // printf("\n\n");
            double**C;
            allocate_matrix_memory(m,p,&C);
            printf("Start of calculations:\n");
            start=clock();
            algorithm_rxsxt_normal(m,n,p,A,B,C,r,s,t,K,P);
            stop=clock();
            printf("Tensors: time %6.3lf\n",(double)(stop-start)/CLOCKS_PER_SEC);

            start=clock();
            ordinary_product(m,n,p,A,B,C);
            stop=clock();
            printf("Ordinary: time %6.3lf\n",(double)(stop-start)/CLOCKS_PER_SEC);

            start=clock();
            algorithm_rxsxt_optim(m,n,p,A,B,C,r,s,t,K,P,L0);
            stop=clock();
            printf("Optimal: time %6.3lf\n",(double)(stop-start)/CLOCKS_PER_SEC);
            printf("\n");

            //print_matrix(m,p,C);

            free_matrix(m,A);
            free_matrix(n,B);
            free_matrix(m,C);
        }
    }
    return 0;
}
int min_integer(int r, int s, int t, int K, int a, int b, int c, double *P[K][3]){
    double x1,x2;
    for(int i=0;i<15;i++){
        x1=f_ord(r,s,t,i);
        x2=f_op(r,s,t,i,K,a,b,c);
        if(x1>x2)
            return i;
    }
    return -100;
}

double f_ord(int r, int s, int t, double l){
    return pow(r,l)*pow(t,l)*(2*pow(s,l)-1);
}

double f_op(int r, int s, int t, double l, int K, int a, int b, int c){
    double sumaord,sumaa,sumab,sumac;
    double k_aux=K;
    sumaord=k_aux*pow(r,l-1)*pow(t,l-1)*(2*pow(s,l-1)-1);
    sumaa=a*pow(r,l-1)*pow(s,l-1);
    sumab=b*pow(s,l-1)*pow(t,l-1);
    sumac=c*pow(r,l-1)*pow(t,l-1);
    return(sumaord+sumaa+sumab+sumac);
}

int additions_tensor_rxs(int r, int s, int K, double *P[K][3]){
    int rxs=0;  //additions of r x s matrix
    int sum;
    //additions r x s:
    for(int k=0;k<K;k++){
        sum=0;
        for(int i=0;i<(r*s);i++){
            if(P[k][0][i]!=0)
                sum++;
        }
        rxs=rxs+(sum-1);  // Number of ones - 1.
    }
    return rxs;
}

int additions_tensor_sxt(int s, int t, int K, double *P[K][3]){
    int sxt=0;  //additions of s x t matrix
    int sum;
    //additions s x t:
    for(int k=0;k<K;k++){
        sum=0;
        for(int i=0;i<(s*t);i++){
            if(P[k][1][i]!=0)
                sum++;
        }
        sxt=sxt+(sum-1);  // Number of ones - 1.
    }
    return sxt;
}

int additions_tensor_rxt(int r, int t, int K, double *P[K][3]){
    int rxt=0;  //additions of r x t matrix
    int sum;
    //additions r x t:
    for(int i=0;i<(r*t);i++){
            sum=0;
        for(int k=0;k<K;k++){
            if(P[k][2][i]!=0)
                sum++;
        }
        rxt=rxt+(sum-1);
    }
    return rxt;
}
                                                    //rows columns
void matrices_combination(int m, int n, double **A, int f, int c, double **A1, double *coef){
    int m1=m/f;
    int n1=n/c;
    matrix_of_zeros(m1,n1,A1);
    for(int i=0;i<f;i++){   //rows of matrices
        for(int j=0;j<c;j++){  //columns of matrices
            if(coef[c*i+j]!=0){
                for(int k=0;k<m1;k++){   //rows of each matrix
                    for(int l=0;l<n1;l++){   //columns of each matrix
                        A1[k][l]=A1[k][l]+A[i*m1+k][j*n1+l]*coef[c*i+j];
                    }
                }
            }
        }
    }
}
                                                    //rows columns
void matrix_recombination(int m, int n, double **C, int f, int c, double **A1, double *coef){
    int m1=m/f;
    int n1=n/c;
    for(int i=0;i<f;i++){   //rows of matrices
        for(int j=0;j<c;j++){  //columns of matrices
            if(coef[c*i+j]!=0){
                for(int k=0;k<m1;k++){   //rows of each matrix
                    for(int l=0;l<n1;l++){   //columns of each matrix
                        C[i*m1+k][j*n1+l]+=A1[k][l]*coef[c*i+j];
                    }
                }
            }
        }
    }
}

void algorithm_rxsxt_normal(int m, int n, int p, double **A, double **B, double **C,int r, int s, int t, int K, double *P[K][3]){
    if(m==1 && m==n && n==p)
        C[0][0]= A[0][0]*B[0][0];
    else{
        matrix_of_zeros(m,p,C);
        double **A1;
        double **B1;
        double **P1;

        allocate_matrix_memory(m/r,n/s,&A1);  //m=r*m1, n=s*n1:  size A1: m1 x n1
        allocate_matrix_memory(n/s,p/t,&B1);  //n=s*n1, p=t*p1   size B1: n1 x p1
        allocate_matrix_memory(m/r,p/t,&P1);  //m=r*m1, p=t*p1   size A1: m1 x p1
        for(int i=0;i<K;i++){
            matrices_combination(m,n,A,r,s,A1,P[i][0]); //Lineal combination of A
            matrices_combination(n,p,B,s,t,B1,P[i][1]);  //Lineal combination of B
            algorithm_rxsxt_normal(m/r,n/s,p/t,A1,B1,P1,r,s,t,K,P);   // A*B in P1
            matrix_recombination(m,p,C,r,t,P1,P[i][2]);  //Pi added in the position of C indicated by P
        }
        free_matrix(m/r,A1);
        free_matrix(n/s,B1);
        free_matrix(m/r,P1);
    }
}


                          //m: rows A, n: columns A and rows B, p: columns B
void algorithm_rxsxt_optim(int m, int n, int p, double **A, double **B, double **C,int r, int s, int t, int K, double *P[K][3],int L0){
    if(m<=pow(r,L0)){
        ordinary_product(m,n,p,A,B,C);
    }
    else{
        matrix_of_zeros(m,p,C);
        double **A1;
        double **B1;
        double **P1;

        allocate_matrix_memory(m/r,n/s,&A1);  //m=r*m1, n=s*n1:  size A1: m1 x n1
        allocate_matrix_memory(n/s,p/t,&B1);  //n=s*n1, p=t*p1   size B1: n1 x p1
        allocate_matrix_memory(m/r,p/t,&P1);  //m=r*m1, p=t*p1   size A1: m1 x p1
        for(int i=0;i<K;i++){
            matrices_combination(m,n,A,r,s,A1,P[i][0]); // Lineal combination of A
            matrices_combination(n,p,B,s,t,B1,P[i][1]);  //Lineal combination of B B
            algorithm_rxsxt_optim(m/r,n/s,p/t,A1,B1,P1,r,s,t,K,P,L0);   //A*B in P1
            matrix_recombination(m,p,C,r,t,P1,P[i][2]);  // Pi added in the C position indicated by P
        }
        free_matrix(m/r,A1);
        free_matrix(n/s,B1);
        free_matrix(m/r,P1);
    }
}


void matrix_of_zeros(int m, int n, double **M){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            M[i][j]=0;
        }
    }
}

void random_matrix(int m, int n, double **M){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            M[i][j]=(double)rand()+(double)rand()/RAND_MAX;
            if(rand()%2==1){
                M[i][j]=-M[i][j];
            }
        }
    }
}

void allocate_matrix_memory(int m, int n, double ***M){
    (*M) = (double **)malloc(m*sizeof(double *));
    for(int i=0;i<m;i++){
        (*M)[i]=(double *)malloc(n*sizeof(double));
    }

    if ((*M) == NULL) {
        printf("Memory allocation failed\n");
    }
}
void free_matrix(int m, double **M) {
    for (int i = 0; i < m; i++) {
        free(M[i]);
    }
    free(M);
}

void print_matrix(int m, int n, double **M){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            printf("%lf   ",M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}


void ordinary_product(int m, int n, int p, double **A, double **B, double **C){
    matrix_of_zeros(m,p,C);
    for(int m_ite=0;m_ite<m;m_ite++){
        for(int p_ite=0;p_ite<p;p_ite++){
            for(int j=0;j<n;j++){
                C[m_ite][p_ite] = C[m_ite][p_ite] + A[m_ite][j]*B[j][p_ite];
            }
        }
    }
}

