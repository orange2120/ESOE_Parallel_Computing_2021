#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DEBUG 0

/**
 * @brief read the matrix from file
 * 
 * @param path path to the matrix binary file
 * @param m rows 
 * @param n cols
 * @param matrix storage of matrix
 */
void readMatrix(char* path, int* m, int* n, double*** matrix) {
    FILE* infileptr = fopen(path, "r");
    if (infileptr == NULL) {
        printf("File %s not found.\n", path);
        exit(1);
    }
    fread(m, sizeof(int), 1, infileptr);
    fread(n, sizeof(int), 1, infileptr);
    *matrix = (double**)malloc(*m * sizeof(void*));
    for (int i = 0; i < *m; ++i) {
        (*matrix)[i] = (double*)malloc(*n * sizeof(double));
        for (int j = 0; j < *n; ++j) {
            fread(&(*matrix)[i][j], sizeof(double), 1, infileptr);
        }
    }
    fclose(infileptr);
}

/**
 * @brief print the matrix
 * 
 * @param mat the matrix
 * @param matM rows
 * @param matN cols
 */
void printMatrix(double** mat, int matM, int matN) {
    for (int i = 0; i < matM; ++i) {
        for (int j = 0; j < matN; ++j) {
            printf("%6.3f ", mat[i][j]);
        }
        printf("\n");
    }
}

/**
 * @brief read the vector
 * 
 * @param path path to the vector file
 * @param n size of the vector
 * @param vector storage of the vector
 */
void readVector(char* path, int* n, double** vector) {
    FILE* infileptr = fopen(path, "r");
    if (infileptr == NULL) {
        printf("File %s not found.\n", path);
        exit(1);
    }
    fread(n, sizeof(int), 1, infileptr);
    *vector = (double*)malloc(*n * sizeof(double));
    for (int i = 0; i < *n; ++i) {
        fread(&(*vector)[i], sizeof(double), 1, infileptr);
    }
    fclose(infileptr);
}

void readOutput(char* path, int* n, double** vector) {
    FILE* infileptr = fopen(path, "r");
    if (infileptr == NULL) {
        printf("File %s not found.\n", path);
        exit(1);
    }
    fscanf(infileptr, "%d", n);
    *vector = (double*)malloc(*n * sizeof(double));
    for (int i = 0; i < *n; ++i) {
        fscanf(infileptr, "%lf", &(*vector)[i]);
    }
    fclose(infileptr);
}

/**
 * @brief print the vector
 * 
 * @param vec storage of the vector
 * @param vecN size of the vector
 */
void printVector(double* vec, int vecN) {
    for (int i = 0; i < vecN; ++i) {
        printf("%6.3f ", vec[i]);
    }
    printf("\n");
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: ./validate <path to matrix> <path to vector> <path to output>\n");
        exit(1);
    }
    // read matrix
    int matM, matN;
    double** mat;
    readMatrix(argv[1], &matM, &matN, &mat);
    printf("[Matrix] rows = %d, cols = %d\n", matM, matN);
    if (DEBUG) printMatrix(mat, matM, matN);
    // read vector
    int vecN;
    double* vec;
    readVector(argv[2], &vecN, &vec);
    printf("[Vector] size = %d\n", vecN);
    if (DEBUG) printVector(vec, vecN);

    // read vector
    int outVecN;
    double* outVec;
    readOutput(argv[3], &outVecN, &outVec);
    printf("[Output Vector] size = %d\n", outVecN);
    if (DEBUG) printVector(outVec, outVecN);

    // validate
    if (matN != vecN) {
        printf("Size not compatible!\n");
        exit(1);
    }

    // calculate matrix * vector
    double* result = (double*)malloc(matM * sizeof(double));
    for (int i = 0; i < matM; ++i) {
        result[i] = 0.0;
        for (int j = 0; j < matN; ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    printf("[Matrix * Vector] size = %d\n", matM);
    if (DEBUG) printVector(result, matM);

    int valid = 1;
    if (matM != outVecN) {
        valid = 0;
    }
    for (int i = 0; i < matM; ++i) {
        // printf("correct = %.8f, yours = %.8f\n", result[i], outVec[i]);
        if (fabs(outVec[i] - result[i]) > 0.0000001f) {
            valid = 0;
        }
    }
    if (!valid) {
        printf("The result is incorrect!\n");

    } else {
        printf("The result is correct!\n");
    }
}