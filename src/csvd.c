

#include <stdio.h>
void print_mat(const double * mat, const int n, const int m){
  printf("\n[");
  for(int i = 0; i < n; i++){
    i? printf(" [") : printf("[");
    for(int j = 0; j < m-1; j ++) printf(" %.10f,", mat[i + n*j]);
    i +1 == n? printf(" %.10f ]]\n", mat[i + n*(m-1)]):printf(" %.10f ],\n", mat[i + n*(m-1)]);
  }
  printf("\n");
}

void mult_row(const double * a, const double *b, const int *nrowA,
              const int *ncolA, const int *ncolB, double * c){
  for(int i = 0; i<*nrowA; i++){
    for(int j = 0; j < *ncolB; j++){
      c[i**ncolB + j] = 0;
      for(int k = 0; k < *ncolA; k++)
        c[i**ncolB + j] += a[i**ncolA + k] * b[j + k**ncolB];
    }
  }
}

void mult_col(const double * a, const double *b, const int *nrowA, const int *ncolA, const int *ncolB, double * c){
  for(int i = 0; i < *nrowA; i++){
    for(int j = 0; j < *ncolB; j++){
      c[j**nrowA + i] = 0;
      for(int k = 0; k < *ncolA; k++)
        c[j**nrowA + i] += a[i + k**nrowA] * b[k + j**ncolA];
    }
  }
}



#include<math.h>


void NORMALIZE(double *x, int n, double * norm){                                      \
  double sum = 0;
  for(int i = 0; i < n; i++)sum += x[i] * x[i];
  double denom = sqrt(sum);
  int d = (isinf(denom) || denom == 0 || isnan(denom));
  if (d) denom = 1;
  for(int i = 0; i < n; i++) x[i]/=denom;
  norm[0] = d? 0: denom;
}


#include<stdlib.h>
void power(const double *a, const int *nrow, double * b){
  double * v = (double*)calloc(*nrow, sizeof(double));
  v[0] = 0.02;
  int ncolB = 1;
  double sd[0];
  int N = *nrow < 30? 10 :(*nrow < 50? 20: 10*log2(*nrow));
  for(int i = 0; i < N; i++){
    mult_col(a, v, nrow, nrow, &ncolB, b);
    for(int j = 0; j < *nrow; j++)v[j] = b[j]/b[0];

  }
  free(v);
  NORMALIZE(b, *nrow, sd);
}






void tcrossprod(const double *A, const int n, const int m, double *B){
  for (int i = 0; i < n; i++)
    for(int j = i; j < n; j++){
      B[i*n + j] = 0;
      for (int k = 0; k < m; k++) B[i*n + j] += A[i + k*n] * A[j + k*n];
      B[i + j*n] = B[i*n +j];
    }
}


void crossprod(const double *A, const int n, const int m, double *B){
  for(int i = 0; i < m; i++)
    for (int j = i; j < m; j++){
      B[i*m + j] = 0;
      for (int k = 0; k < n; k++) B[i*m + j] +=  A[i*n + k] * A[j*n + k];
      B[i+j*m] = B[i*m+j];
    }
}

void mult_col_sub(const double * a, const double *b, const int *nrowA,
                  const int *ncolA, const int *ncolB, double * c, double S){
  for(int i = 0; i < *nrowA; i++){
    for(int j = 0; j < *ncolB; j++){
      for(int k = 0; k < *ncolA; k++)
        c[j**nrowA + i] -= S*a[i + k**nrowA] * b[k + j**ncolA];

    }
  }
}


void SV(double *A, const int *nrowA, const int *ncolA,  double*S,double *eigenU, double*eigenV){
  int transposed = (*nrowA < *ncolA);
  void(*cross)(const double*,const int, const int, double*);
  void(*mult)(const double*,const double*, const int*,const int*,const int*, double*);
  if(transposed) {
    cross = tcrossprod;
    mult = mult_row;
  }
  else {
    cross = crossprod;
    mult = mult_col;
  }

  int NCOL = transposed? *nrowA : *ncolA;
  int NROW = transposed? *ncolA : *nrowA;
  double *D = (double*)calloc(NCOL*NCOL, sizeof(double));
  int ncolB = 1;
  for(int i = 0; i < NCOL; i++){
    double * eigenVptr = eigenV + i * NCOL;
    double * eigenUptr = eigenU + i * NROW;
    cross(A, *nrowA, *ncolA,  D);
    power(D, &NCOL, eigenVptr);
    mult(A, eigenVptr, &NROW, &NCOL, &ncolB, eigenUptr);
    NORMALIZE(eigenUptr, NROW, S+i);
    if(transposed)
      mult_col_sub(eigenVptr,  eigenUptr, nrowA, &ncolB, ncolA, A, S[i]);
    else
      mult_col_sub(eigenUptr,  eigenVptr, nrowA, &ncolB, ncolA, A, S[i]);

  }
  free(D);
}













int main(){
  double a[15]={1,2,3,4,5,6,7,8,9, 10,11,12,13,14,15};
  //double d[50];
  int n = 3;

  int mm = 5;
  //double b[15] = {1, 4, 7, 10, 13, 2, 5, 8, 11, 14, 3, 6, 9, 12, 15};

  double S[3];
  double U[15];
  double V[9];
  SV(a, &n, &mm, S, U, V);
  print_mat(S, 1, 3);
  print_mat(U, 5, 3);
  print_mat(V, 3, 3);
  //SV(b, &mm, &n);
}
