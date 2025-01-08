#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void random_matrix(int m, int n,FILE *matrix){
    double aux;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            aux=(double)rand()+(double)rand()/RAND_MAX;
            if(rand()%2==1)
                fprintf(matrix,"%lf ",-aux);
            else
                fprintf(matrix,"%lf ",aux);
        }
        fprintf(matrix,"\n");
    }
}

int main(){

    int m, n, p;
    printf("m: ");
    scanf("%d",&m);
    printf("n: ");
    scanf("%d",&n);
    printf("p: ");
    scanf("%d",&p);

    double **A = (double **)malloc(m * sizeof(double *));
    if (A == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    for (int i = 0; i < m; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
        if (A[i] == NULL) {
            printf("Memory allocation failed\n");
            for (int j = 0; j < i; j++) {
                free(A[j]);
            }
            free(A);
            return 1;
        }
    }

    char fitxer_A[20];
    printf("File name for first matrix: ");
    scanf("%s",fitxer_A);
    FILE *matrixA = fopen(fitxer_A, "w");
    if (matrixA == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    random_matrix(m,n,matrixA);
    fclose(matrixA);
    for (int i = 0; i < m; i++) {
        free(A[i]);
    }
    free(A);

    printf("A done\n");

    double **B = (double **)malloc(n * sizeof(double *));
    if (B == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    for (int i = 0; i < n; i++) {
        B[i] = (double *)malloc(p * sizeof(double));
        if (B[i] == NULL) {
            printf("Memory allocation failed\n");
            for (int j = 0; j < i; j++) {
                free(B[j]);
            }
            free(B);
            return 1;
        }
    }
    char fitxer_B[20];
    printf("File name for second matrix: ");
    scanf("%s",fitxer_B);
    FILE *matrixB = fopen(fitxer_B, "w");
    if (matrixB == NULL) {
        printf("No se puede abrir el archivo.\n");
        return 1;
    }

    random_matrix(n,p,matrixB);

    fclose(matrixB);
     for (int i = 0; i < n; i++) {
        free(B[i]);
    }
    free(B);
    return 0;
}
