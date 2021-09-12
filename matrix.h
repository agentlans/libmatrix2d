#ifndef _MATRIX
#define _MATRIX

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. 
 *
 * Copyright 2021 Alan Tseng */

// Matrix in row-major order
typedef struct matrix_t {
    double* data;
    int rows;
    int cols;
} matrix_t;

// To create and destroy matrix
matrix_t* matrix_new(int rows, int cols);
void matrix_free(matrix_t* m);

// Getting and setting individual elements of matrix
double matrix_get(matrix_t* m, int row, int col);
void matrix_set(matrix_t* m, int row, int col, double x);
// Returns pointer to the matrix element
double* matrix_ref(matrix_t* m, int row, int col);

// Input and output from files
matrix_t* matrix_read(const char* filename);
int matrix_write(const char* filename, matrix_t* m);

// Stacking matrices by row or column
matrix_t* matrix_bind_rows(matrix_t** matrices, int num_matrices);
matrix_t* matrix_bind_cols(matrix_t** matrices, int num_matrices);

// Extract a single row or column
void matrix_get_row(double* out, matrix_t* m, int row);
void matrix_get_col(double* out, matrix_t* m, int col);

// Overwrite a single row or column
void matrix_set_row(matrix_t* m, int row, double* replacement);
void matrix_set_col(matrix_t* m, int col, double* replacement);

// Return a submatrix containing multiple rows or columns of the original matrix
matrix_t* matrix_subset_rows(matrix_t* m, int* row_id, int num_rows);
matrix_t* matrix_subset_cols(matrix_t* m, int* col_id, int num_cols);
matrix_t* matrix_subset(matrix_t* m, int* row_id, int num_rows, int* col_id, int num_cols);

void matrix_print(matrix_t* m);
matrix_t* matrix_transpose(matrix_t* m);

// Fill matrix with a number or consecutive numbers
void matrix_fill(matrix_t* m, double x);
void matrix_fill_consecutive(matrix_t* m);

// Matrix diagonals
void matrix_get_diagonal(double* out, matrix_t* m);
void matrix_set_diagonal(matrix_t* m, double* out);

// Copies a matrix to dest[row,col], overwriting previous values
void matrix_paste(matrix_t* dest, int row, int col, matrix_t* src);

#endif
