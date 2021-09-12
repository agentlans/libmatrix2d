#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "matrix.h"

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. 
 *
 * Copyright 2021 Alan Tseng */

// Returns pointer to a rows x cols matrix
matrix_t* matrix_new(int rows, int cols) {
    long n = rows * cols;
    double* data = calloc(n, sizeof(double));
    matrix_t* m = malloc(sizeof(matrix_t));
    if (!data || !m) {
        free(data);
        free(m);
        return NULL;
    }
    *m = (matrix_t) {data, rows, cols};
    return m;
}

void matrix_free(matrix_t* m) {
    if (m) free(m->data);
    free(m);
}

// Returns pointer to m[row,col]
double* matrix_ref(matrix_t* m, int row, int col) {
    return m->data + row*m->cols + col;
}

// Returns m[row,col]
double matrix_get(matrix_t* m, int row, int col) {
    return *matrix_ref(m, row, col);
}

// Sets m[row,col] = x
void matrix_set(matrix_t* m, int row, int col, double x) {
    *(matrix_ref(m, row, col)) = x;
}

// Copies n doubles from src to dest
void copy_doubles(double* dest, double* src, int n) {
    memmove(dest, src, n*sizeof(double));
}

// Copies row of matrix to an array
void matrix_get_row(double* out, matrix_t* m, int row) {
    copy_doubles(out, matrix_ref(m, row, 0), m->cols);
}

// Copies array into the row of a matrix
void matrix_set_row(matrix_t* m, int row, double* replacement) {
    copy_doubles(matrix_ref(m, row, 0), replacement, m->cols);
}

// Copies column of matrix to an array
void matrix_get_col(double* out, matrix_t* m, int col) {
    int r = m->rows;
    for (int i = 0; i < r; ++i) {
        out[i] = matrix_get(m, i, col);
    }
}

// Copies array into the column of a matrix
void matrix_set_col(matrix_t* m, int col, double* replacement) {
    int r = m->rows;
    for (int i = 0; i < r; ++i) {
        matrix_set(m, i, col, replacement[i]);
    }
}

// Returns a matrix containing rows of the original matrix
matrix_t* matrix_subset_rows(matrix_t* m, int* row_id, int num_rows) {
    matrix_t* out = matrix_new(num_rows, m->cols);
    if (!out) return NULL;
    double tmp[m->cols]; // Store the row temporarily
    for (int i = 0; i < num_rows; ++i) {
        matrix_get_row(tmp, m, row_id[i]);
        matrix_set_row(out, i, tmp);
    }
    return out;
}

// Returns a matrix containing columns of the original matrix
matrix_t* matrix_subset_cols(matrix_t* m, int* col_id, int num_cols) {
    matrix_t* out = matrix_new(m->rows, num_cols);
    if (!out) return NULL;
    double tmp[m->rows]; // Store the column temporarily
    for (int i = 0; i < num_cols; ++i) {
        matrix_get_col(tmp, m, col_id[i]);
        matrix_set_col(out, i, tmp);
    }
    return out;
}

// Returns a matrix containing a subset of the rows and columns of the original matrix
matrix_t* matrix_subset(matrix_t* m, int* row_id, int num_rows, int* col_id, int num_cols) {
    matrix_t* out = matrix_new(num_rows, num_cols);
    if (!out) return NULL;
    // Go over the indices
    for (int i = 0; i < num_rows; ++i) {
        int row = row_id[i];
        for (int j = 0; j < num_cols; ++j) {
            int col = col_id[j];
            // Copy original element to the output matrix
            double x = matrix_get(m, row, col);
            matrix_set(out, i, j, x);
        }
    }
    return out;
}

// Pastes the src matrix at dest[row,col]
void matrix_paste(matrix_t* dest, int row, int col, matrix_t* src) {
    for (int i = 0; i < src->rows; ++i) {
        int dest_row = row + i;
        for (int j = 0; j < src->cols; ++j) {
            int dest_col = col + j;
            double x = matrix_get(src, i, j);
            matrix_set(dest, dest_row, dest_col, x);
        }
    }
}

matrix_t* matrix_transpose(matrix_t* m) {
    matrix_t* out = matrix_new(m->cols, m->rows);
    if (!out) return NULL;
    for (int i = 0; i < m->rows; ++i) {
        for (int j = 0; j < m->cols; ++j) {
            double x = matrix_get(m, i, j);
            matrix_set(out, j, i, x);
        }
    }
    return out;
}

// Fills matrix with the value x
void matrix_fill(matrix_t* m, double x) {
    long n = m->rows * m->cols;
    for (long i = 0; i < n; ++i) {
        m->data[i] = x;
    }
}

// Fills matrix with consecutive integers starting from 0
void matrix_fill_consecutive(matrix_t* m) {
    long n = m->rows * m->cols;
    for (long i = 0; i < n; ++i) {
        m->data[i] = i;
    }
}

int min(int a, int b) {
    if (a < b) {
        return a;
    } else {
        return b;
    }
}

int smallest_dim(matrix_t* m) {
    return min(m->rows, m->cols);
}

void matrix_get_diagonal(double* out, matrix_t* m) {
    int n = smallest_dim(m);
    for (int i = 0; i < n; ++i) {
        out[i] = matrix_get(m, i, i);
    }
}

void matrix_set_diagonal(matrix_t* m, double* out) {
    int n = smallest_dim(m);
    for (int i = 0; i < n; ++i) {
        matrix_set(m, i, i, out[i]);
    }
}

// Stack matrices by rows
matrix_t* matrix_bind_rows(matrix_t** matrices, int num_matrices) {
    if (num_matrices == 0) return NULL;
    // Figure out the size of the stacked matrix
    int total_rows = 0;
    int cols = -1;
    for (int i = 0; i < num_matrices; ++i) {
        matrix_t* mat = matrices[i];
        if (cols == -1) {
            cols = mat->cols;
        } else {
            // The matrices must have the same number of columns
            assert(cols == mat->cols);
        }
        total_rows += mat->rows;
    }
    // The stacked matrix
    matrix_t* stacked = matrix_new(total_rows, cols);
    if (!stacked) return NULL;
    // Stack rows
    int start = 0;
    for (int i = 0; i < num_matrices; ++i) {
        matrix_t* m = matrices[i];
        matrix_paste(stacked, start, 0, m);
        start += m->rows;
    }
    return stacked;
}

// Stack matrices by columns
matrix_t* matrix_bind_cols(matrix_t** matrices, int num_matrices) {
    if (num_matrices == 0) return NULL;
    // Figure out the size of the stacked matrix
    int rows = -1;
    int total_cols = 0;
    for (int i = 0; i < num_matrices; ++i) {
        matrix_t* mat = matrices[i];
        if (rows == -1) {
            rows = mat->rows;
        } else {
            // The matrices must have the same number of rows
            assert(rows == mat->rows);
        }
        total_cols += mat->cols;
    }
    // The stacked matrix
    matrix_t* stacked = matrix_new(rows, total_cols);
    if (!stacked) return NULL;
    // Stack columns
    int start = 0;
    for (int i = 0; i < num_matrices; ++i) {
        matrix_t* m = matrices[i];
        matrix_paste(stacked, 0, start, m);
        start += m->cols;
    }
    return stacked;
}

// Prints matrix to stdout
void matrix_print(matrix_t* m) {
    for (int i = 0; i < m->rows; ++i) {
        for (int j = 0; j < m->cols; ++j) {
            printf("%lf", matrix_get(m, i, j));
            if (j != m->cols - 1) {
                printf(" ");
            }
        }
        printf("\n");
    }
}

// Writes matrix to file
int matrix_write(const char* filename, matrix_t* m) {
    FILE* f = fopen(filename, "w");
    if (!f) return 1; // Can't open file
    // Write the dimensions of matrix
    fprintf(f, "%d %d\n", m->rows, m->cols);
    // Write the contents of the matrix
    for (int i = 0; i < m->rows; ++i) {
        for (int j = 0; j < m->cols; ++j) {
            fprintf(f, "%lf", matrix_get(m, i, j));
            if (j != m->cols - 1) {
                fprintf(f, " ");
            }
        }
        fprintf(f, "\n");
    }
    fclose(f);
    return 0;
}

// Reads matrix from file
matrix_t* matrix_read(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) return NULL; // Can't open file
    // Read dimensions of matrix
    int rows, cols;
    int res = fscanf(f, "%d %d\n", &rows, &cols);
    assert(res == 2);
    // Set space for matrix
    matrix_t* m = matrix_new(rows, cols);
    if (!m) {
        fclose(f);
        return NULL;
    }
    // Read the file and save into matrix
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double x;
            int res = fscanf(f, "%lf", &x);
            assert(res == 1);
            matrix_set(m, i, j, x);
        }
    }
    fclose(f);
    return m;
}

