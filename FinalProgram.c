#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <malloc.h>
int calculate_l(int m_ini,int n_ini,int p_ini,int Ti[3]);
int max3(int l1,int l2,int l3);
int find_l(double arr);
void random_matrix_to_tensor_matrix(int m_ini,int n_ini,int p_ini,int T[93][3], int rst[4]);

void matrix_of_zeros(int m, int n, double **M);
void allocate_matrix_memory(int m, int n, double ***M);
void matrix_redimension(int m_ini, int n_ini, int m, int n, double ***M);
void ordinary_product(int m, int n, int p, double **A, double **B, double **C);
void submatrices_combination(int m, int n, double **A, int f, int c, double **A1, double *coef);
void submatrix_recombination(int m, int n, double **C, int f, int c, double **A1, double *coef);
void algorithm_rxsxt(int m, int n, int p, double **A, double **B, double **C,int r, int s, int t, int K, double *P[K][3],int L0);
void free_matrix(int m, double **M);

int additions_tensor_rxs(int r, int s, int K, double *P[K][3]);
int additions_tensor_sxt(int s, int t, int K, double *P[K][3]);
int additions_tensor_rxt(int r, int t, int K, double *P[K][3]);

double f_ord(int r, int s, int t, double l);
double f_op(int r, int s, int t, double l, int K, int a, int b, int c);

int min_integer_optimal(int r, int s, int t, int K, int a, int b, int c, double *P[K][3]);
int minimZero(int m_ini,int n_ini,int p_ini,int L0_op[93],int aux[5]);
void optimal_threshold(char *tensor_names[93],int L0_op[93]);

char *tensor_names[93] = {"tensor_2x2x2.txt",
                            "tensor_2x2x3.txt",
                            "tensor_2x2x4.txt",
                            "tensor_2x2x5.txt",
                            "tensor_2x2x6.txt",
                            "tensor_2x2x7.txt",
                            "tensor_2x2x8.txt",
                            "tensor_2x3x3.txt",
                            "tensor_2x3x4.txt",
                            "tensor_2x3x5.txt",
                            "tensor_2x4x4.txt",
                            "tensor_2x4x5.txt",
                            "tensor_2x5x5.txt",
                            "tensor_3x3x3.txt",
                            "tensor_3x3x4.txt",
                            "tensor_3x3x5.txt",
                            "tensor_3x4x11.txt",
                            "tensor_3x4x4.txt",
                            "tensor_3x4x5.txt",
                            "tensor_3x5x5.txt",
                            "tensor_3x5x9.txt",
                            "tensor_3x9x11.txt",
                            "tensor_4x4x4.txt",
                            "tensor_4x4x5.txt",
                            "tensor_4x5x10.txt",
                            "tensor_4x5x11.txt",
                            "tensor_4x5x5.txt",
                            "tensor_4x5x9.txt",
                            "tensor_4x9x10.txt",
                            "tensor_4x9x11.txt",
                            "tensor_4x11x11.txt",
                            "tensor_4x11x12.txt",
                            "tensor_5x5x5.txt",
                            "tensor_5x5x7.txt",
                            "tensor_5x7x10.txt",
                            "tensor_5x7x11.txt",
                            "tensor_5x7x9.txt",
                            "tensor_5x8x10.txt",
                            "tensor_5x8x11.txt",
                            "tensor_5x8x9.txt",
                            "tensor_5x9x10.txt",
                            "tensor_5x9x11.txt",
                            "tensor_5x9x12.txt",
                            "tensor_5x9x9.txt",
                            "tensor_6x7x10.txt",
                            "tensor_6x7x11.txt",
                            "tensor_6x7x9.txt",
                            "tensor_6x8x10.txt",
                            "tensor_6x8x11.txt",
                            "tensor_6x9x10.txt",
                            "tensor_6x9x11.txt",
                            "tensor_6x9x9.txt",
                            "tensor_7x7x10.txt",
                            "tensor_7x7x11.txt",
                            "tensor_7x7x9.txt",
                            "tensor_7x8x10.txt",
                            "tensor_7x8x11.txt",
                            "tensor_7x8x12.txt",
                            "tensor_7x8x9.txt",
                            "tensor_7x9x10.txt",
                            "tensor_7x9x11.txt",
                            "tensor_7x9x12.txt",
                            "tensor_7x9x9.txt",
                            "tensor_7x10x10.txt",
                            "tensor_7x10x11.txt",
                            "tensor_7x11x11.txt",
                            "tensor_8x8x10.txt",
                            "tensor_8x8x11.txt",
                            "tensor_8x9x10.txt",
                            "tensor_8x9x11.txt",
                            "tensor_8x9x12.txt",
                            "tensor_8x10x10.txt",
                            "tensor_8x10x11.txt",
                            "tensor_8x10x12.txt",
                            "tensor_8x11x11.txt",
                            "tensor_8x11x12.txt",
                            "tensor_9x9x10.txt",
                            "tensor_9x9x11.txt",
                            "tensor_9x9x9.txt",
                            "tensor_9x10x10.txt",
                            "tensor_9x10x11.txt",
                            "tensor_9x10x12.txt",
                            "tensor_9x11x11.txt",
                            "tensor_9x11x12.txt",
                            "tensor_10x10x10.txt",
                            "tensor_10x10x11.txt",
                            "tensor_10x10x12.txt",
                            "tensor_10x11x11.txt",
                            "tensor_10x11x12.txt",
                            "tensor_10x12x12.txt",
                            "tensor_11x11x11.txt",
                            "tensor_11x11x12.txt",
                            "tensor_11x12x12.txt"};


int main(){

    ///READING OF MATRIX A
    FILE *A_txt;

    char A_file[20];

    printf("Matrix A: ");
    scanf("%s",A_file);
    A_txt = fopen(A_file, "r");

    if(A_txt == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    int m_ini=0;
    int n_ini=0;
    char c;
    while((c = fgetc(A_txt)) != EOF) {
        if(c == '\n')
            m_ini++; // Counts each \n
    }
    rewind(A_txt);
    while((c = fgetc(A_txt)) != '\n'){
        if(c == ' ') // Counts the black spaces
            n_ini++;
    }
    rewind(A_txt);

    double **A;
    allocate_matrix_memory(m_ini,n_ini,&A);

    for (int i=0;i<m_ini;i++) {
        for (int j=0;j<n_ini;j++)
            fscanf(A_txt, "%lf", &A[i][j]);
    }
    fclose(A_txt);

    ///READING MATRIX B
    FILE *B_txt;

    char B_file[20];

    printf("Matrix B: ");
    scanf("%s",B_file);
    B_txt = fopen(B_file, "r");

    if(B_txt == NULL) {
        printf("Error opening file\n");
        return 1;
    }
    int n_ini_B=0;
    int p_ini=0;
    while((c = fgetc(B_txt)) != EOF) {
        if(c == '\n')
            n_ini_B++; // counts each \n
    }
    if(n_ini != n_ini_B){
        printf("The number of columns of A does not match the number of rows if B\n");
        return 1;
    }
    rewind(B_txt);
    while((c = fgetc(B_txt)) != '\n'){
        if(c == ' ') //Conuts the blank spaces
            p_ini++;
    }
    rewind(B_txt);

    double **B;
    allocate_matrix_memory(n_ini,p_ini,&B);

    for (int i=0;i<n_ini;i++) {
        for (int j=0;j<p_ini;j++)
            fscanf(B_txt, "%lf", &B[i][j]);
    }
    fclose(B_txt);

    printf("Initial matrix multiplication of size %dx%dx%d\n",m_ini,n_ini,p_ini);

    clock_t start, stop;

    ///OPTIMAL THRESHOLD OF EACH TENSOR
    int L0_op[93]={5,4,4,5,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,4,3,3,3,3,3,3,
                   3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                   3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                   3,3,3,3,3,3,3,3,3};
    ///ZEROS OF EACH TENSOR
    // Choses the tensor with the minimal number of zeros
    int aux[5];
    minimZero(m_ini,n_ini,p_ini,L0_op,aux);
    int r,s,t,l,L0;
    r=aux[0];
    s=aux[1];
    t=aux[2];
    l=aux[3];
    L0=aux[4];

    if(l<L0){ //if the power needed l is smaller that the threshold: ordinary multiplication
        double**C;
        allocate_matrix_memory(m_ini,p_ini,&C);
        start=clock();
        ordinary_product(m_ini,n_ini,p_ini,A,B,C);
        stop=clock();
        printf("The best option is the ordinary product. Time:%6.3lf\n",(double)(stop - start) / CLOCKS_PER_SEC);
    }
    else{ // If l>threshold, tensor product
        int m=pow(r,l);
        int n=pow(s,l);
        int p=pow(t,l);

        ///ORDINARY PRODUCT
        double**C_ord;
        allocate_matrix_memory(m_ini,p_ini,&C_ord);

        start = clock();
        ordinary_product(m_ini,n_ini,p_ini,A,B,C_ord);
        stop = clock();
        printf("Ordinary multiplcation time: %6.3lf\n",(double)(stop - start) / CLOCKS_PER_SEC);
        free_matrix(m_ini,C_ord);

        ///ADD ZEROS TO THE NEW MATRIX

       matrix_redimension(m_ini,n_ini,m,n,&A);
       matrix_redimension(n_ini,p_ini,n,p,&B);

        double**C;
        allocate_matrix_memory(m,p,&C);

        /// TENSOR READING
        FILE *tensor;

        char tensor_file[20];
        sprintf(tensor_file, "tensor_%dx%dx%d.txt", r, s, t);
        tensor = fopen(tensor_file, "r");

        if(tensor == NULL) {
            printf("Error opening file\n");
            return 1;
        }
        int K;
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

        /// COMPUTE THE TENSOR MULTIPLICATION
        start = clock();
        algorithm_rxsxt(m,n,p,A,B,C,r,s,t,K,P,L0);
        stop = clock();
        printf("Tensor multiplication time: %6.3lf\n",(double)(stop - start) / CLOCKS_PER_SEC);

        free_matrix(m,A);
        free_matrix(n,B);
        free_matrix(m,C);

        for(int i=0;i<K;i++){
            free(P[i][0]);
            free(P[i][1]);
            free(P[i][2]);
        }
    }

    return 0;
}
int minimZero(int m_ini,int n_ini,int p_ini,int L0_op[93],int aux[5]){
    //T[i][1,2,3] are the corresponding r s i t fixated of the 93 tensors
    int T[93][3];
    for(int i=0;i<93;i++)
        sscanf(tensor_names[i],"tensor_%dx%dx%d.txt", &T[i][0], &T[i][1], &T[i][2]);

    long int Z[93];
    int l[93];
    int z1,z2,z3;
    ///CALCULATION OF l AND THE NUMBER OF ADDED ZEROS NEEDED
    for(int i=0;i<93;i++){
        l[i]=calculate_l(m_ini,n_ini,p_ini,T[i]);

        int m_aux=pow(T[i][0],l[i]);
        int n_aux=pow(T[i][1],l[i]);
        int p_aux=pow(T[i][2],l[i]);
        if(m_aux<0 || n_aux<0 || p_aux<0){ //IF m n p <0 the power is too big
            Z[i]=-1;
        }
        else{
            z1=pow(T[i][0],l[i])-m_ini;
            z2=pow(T[i][1],l[i])-n_ini;
            z3=pow(T[i][2],l[i])-p_ini;
            //Each row of zeros counted is multiplied by the number of columns that need to be filled of zeros and vice versa
                 //rows A   columns A      rows B   columns B      rows C   columns C
            Z[i]=z1*n_ini + z2*(m_ini+z1)+ z2*p_ini + z3*(n_ini+z2) + z1*p_ini+ z3*(m_ini+z1);
        }
    }
    ///FIND THE POSITION OF THE TENSOR WITH THE MINIMUM NUMBER OF NEEDED ZEROS
    int row=0;
    while(Z[row]<0){
        row=row+1;
    }

    for(int i=0;i<93;i++){
        if(Z[i]>=0){
            if(Z[i]<Z[row]){
                row=i;
            }
        }
    }
    printf("The optimal tensor to choose for %d,%d,%d is %d,%d,%d and l=%d\n",m_ini,n_ini,p_ini,T[row][0],T[row][1],T[row][2],l[row]);
    aux[0]=T[row][0];
    aux[1]=T[row][1];
    aux[2]=T[row][2];
    aux[3]=l[row];
    aux[4]=L0_op[row];
}

int calculate_l(int m_ini,int n_ini,int p_ini,int Ti[3]){
    double arr1,arr2,arr3;

    arr1=(log(m_ini)/log(Ti[0]));
    arr2=(log(n_ini)/log(Ti[1]));
    arr3=(log(p_ini)/log(Ti[2]));

    int l,l1,l2,l3;
    //With the power arr in't enough, because it's smaller than the original number. We need the following power
    l1=find_l(arr1);
    l2=find_l(arr2);
    l3=find_l(arr3);

    //We grab the bigger l to find the columns and rows of zeros to add
    //The bigger because is the minimum for r,s or t, and we need the matrix to be bigger than all of the original
    return max3(l1,l2,l3);
}

int max3(int l1,int l2,int l3){
    if(l1>=l2)
        if(l1>=l3)
            return l1;
        else
            return l3;
    else
        if(l2>=l3)
            return l2;
        else
            return l3;
}

int find_l(double arr){
    if((arr-(int)arr)<0.00001) //Is an exact power
        return (int)arr;
    else
        return (int)(arr+1);

}



                                                     //rows columns
void submatrices_combination(int m, int n, double **A, int f, int c, double **A1, double *coef){
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
void submatrix_recombination(int m, int n, double **C, int f, int c, double **A1, double *coef){
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



 //m: rows A, n: columns A i rows B, p: columns B
void algorithm_rxsxt(int m, int n, int p, double **A, double **B, double **C,int r, int s, int t, int K, double *P[K][3],int L0){
    if(m<=pow(r,L0)){
        ordinary_product(m,n,p,A,B,C);
    }
    else{
        matrix_of_zeros(m,p,C);
        double **A1;
        double **B1;
        double **P1;

        allocate_matrix_memory(m/r,n/s,&A1);  //m=r*m1, n=s*n1:  sizes A1: m1 x n1
        allocate_matrix_memory(n/s,p/t,&B1);  //n=s*n1, p=t*p1   sizes B1: n1 x p1
        allocate_matrix_memory(m/r,p/t,&P1);  //m=r*m1, p=t*p1   sizes A1: m1 x p1
        for(int i=0;i<K;i++){
            submatrices_combination(m,n,A,r,s,A1,P[i][0]); //lineal combination of A saved in A1
            submatrices_combination(n,p,B,s,t,B1,P[i][1]);  //lineal combination of B saved in B1
            algorithm_rxsxt(m/r,n/s,p/t,A1,B1,P1,r,s,t,K,P,L0);   //multiplication: A1*B1
            submatrix_recombination(m,p,C,r,t,P1,P[i][2]);  //Add P1 ("P_i") in the corresponding C position
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


void allocate_matrix_memory(int m, int n, double ***M){
    (*M) = (double **)malloc(m*sizeof(double *));
    if (*M == NULL) {
        printf("Memory reallocation failed\n");
        exit(1);
    }
    for(int i=0;i<m;i++){
        (*M)[i]=(double *)malloc(n*sizeof(double));
    }

    if ((*M) == NULL) {
        printf("Memory allocation failed\n");
    }
}

void matrix_redimension(int m_ini, int n_ini, int m, int n, double ***M) {
    // Redimensionar les files
    *M = (double **)realloc(*M, m * sizeof(double *));
    if (*M == NULL) {
        printf("Memory reallocation failed\n");
        exit(1);
    }

    // Redimensionar cada fila
    for (int i = 0; i < m; i++) {
        if (i < m_ini) {
            (*M)[i] = (double *)realloc((*M)[i], n*sizeof(double));
        } else {
            (*M)[i] = (double *)malloc(n*sizeof(double));
        }
        if ((*M)[i] == NULL) {
            printf("Memory reallocation failed\n");
            exit(1);
        }
    }

    // Inicializar los nuevos elementos a 0
    for (int i = 0; i < m; i++) {
        for (int j = n_ini; j < n; j++) {
            (*M)[i][j] = 0.0;
        }
    }
    for (int i = m_ini; i < m; i++) {
        for (int j = 0; j < n; j++) {  //podria ser fins n_ini, pero no ens arrisquem...
            (*M)[i][j] = 0.0;
        }
    }
}
                 //free
void free_matrix(int m, double **M) {
    for (int i = 0; i < m; i++) {
        free(M[i]);
    }
    free(M);
}

void imprimir_matriu(int m, int n, double **M){
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


///threshold FUNCTIONS


void optimal_threshold(char *tensor_names[93],int L0_op[93]){

    int cnt = 0;
    while(cnt<93) {
        // Open the file for reading
        FILE *tensor = fopen(tensor_names[cnt], "r");
        if(tensor == NULL) {
            printf("Error opening file\n");
            //return 1;
            continue;
        }
        int r,s,t;
        sscanf(tensor_names[cnt],"tensor_%dx%dx%d.txt", &r, &s, &t);
        int K;
        fscanf(tensor,"%d",&K);
        double *P[K][3];
        for(int i=0;i<K;i++){
            P[i][0]=(double*)malloc((r*s)*sizeof(double));
            P[i][1]=(double*)malloc((s*t)*sizeof(double));
            P[i][2]=(double*)malloc((r*t)*sizeof(double));
            if(P[i][0] == NULL || P[i][1] == NULL || P[i][2] == NULL){
                printf("Memory allocation failed\n");
            }
        }

        for (int i=0;i<K;i++) {
            for (int k=0;k<r*s;k++)
                fscanf(tensor,"%lf", &P[i][0][k]);
            for (int k=0;k<s*t;k++)
                fscanf(tensor,"%lf", &P[i][1][k]);
            for (int k=0;k<r*t;k++)
                fscanf(tensor,"%lf", &P[i][2][k]);
        }
        fclose(tensor);

        int a,b,c;

        a=additions_tensor_rxs(r,s,K,P);
        b=additions_tensor_sxt(s,t,K,P);
        c=additions_tensor_rxt(r,t,K,P);
        L0_op[cnt]=min_integer_optimal(r,s,t,K,a,b,c,P);

        // Free memory
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < 3; j++) {
                free(P[i][j]);
            }
        }
        cnt++;
    }
}

int min_integer_optimal(int r, int s, int t, int K, int a, int b, int c, double *P[K][3]){
    double x1,x2;
    for(int i=0;i<15;i++){
        x1=f_ord(r,s,t,(double)i);
        x2=f_op(r,s,t,(double)i,K,a,b,c);
        //printf("%d ord:%lf >? optim:%lf\n\n",i,x1,x2);
        if(x1>x2){
            return i;
        }
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
    int rxs=0;  //additions de matrius r x s
    int sum;
    //additions r x s:
    for(int k=0;k<K;k++){
        sum=0;
        for(int i=0;i<(r*s);i++){
            if(P[k][0][i]!=0)
                sum++;
        }
        rxs=rxs+(sum-1);  // A cada matriu hi ha el valor de uns menys 1.
    }
    return rxs;
}

int additions_tensor_sxt(int s, int t, int K, double *P[K][3]){
    int sxt=0;  //additions de matrius s x t
    int sum;
    //additions s x t:
    for(int k=0;k<K;k++){
        sum=0;
        for(int i=0;i<(s*t);i++){
            if(P[k][1][i]!=0)
                sum++;
        }
        sxt=sxt+(sum-1);  // A cada matriu hi ha el valor de uns menys 1.
    }
    return sxt;
}

int additions_tensor_rxt(int r, int t, int K, double *P[K][3]){
    int rxt=0;  //additions de matrius r x t
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


