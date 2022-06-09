/*
 *  Created by Hien Vu on 04/21/2021.
 *  Based on MatrixMath.h for Arduino
 */

#include "matrix.h"
#include <math.h>

#define NR_END 1

/* Matrix Printing Routine
// Uses tabs to separate numbers under assumption printed float_prec width won't cause problems
void MatrixPrint(float_prec* A, int m, int n, String label)
{
	// A = input matrix (m x n)
	int i, j;
	//Serial.println();
	//Serial.println(label);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			//Serial.print(A[n * i + j],10);
			//Serial.print("\t");
		}
		//Serial.println();
	}
}
*/

void matrix_error_handle (void){
	//SetPin(Red_LD_GPIO_Port, Red_LD_Pin, SetIO);
	while(1);
}

void matrix_error_handle_inv (int led_on){
	//if (led_on) SetPin(Red_LD_GPIO_Port, Red_LD_Pin, SetIO);
	//else SetPin(Red_LD_GPIO_Port, Red_LD_Pin, ResetIO);
}

// Copy A <= B
void MatrixCopy(float_prec* A, float_prec* B, int row_A, int col_A){
	int col, row;
	for(row = 0; row < row_A; row++){
		for (col = 0; col < col_A; col++){
			A[col_A * row + col] = B[col_A * row + col];
		}
	}
}

// Edited: scan by row first
// A <= Data
void MatrixWrite(float_prec* A, int row_A, int col_A, float_prec data[]){
	int col, row;
	for(row = 0; row < row_A; row++){
		for (col = 0; col < col_A; col++){
			A[col_A * row + col] = data[col_A * row + col];
		}
	}
}

/* Original: scan by column first
// A = input matrix (n x m)
void MatrixWrite(float_prec* A, int row_A, int col_A, float_prec data[]){
  int col, row;
  for (col = 0; col < col_A; col++)
    for(row = 0; row < row_A; row++)
    {
      A[row_A * col + row] = data[row_A * col + row];
    }
}
*/

// A = input matrix (n x n)
void MatrixIdentity(float_prec* A, int dimension){
	int row, col;
	for(row = 0; row < dimension; row++){
  		for (col = 0; col < dimension; col++){
  			if ( row == col )
  				A[dimension * row + col] = 1;
  			else
  				A[dimension * row + col] = 0;
  		}
  	}
}

// A = input matrix (m x n)
// m = number of rows in A
// p = number of columns in A
void MatrixZeros(float_prec* A, int n, int m){
  int i, j;
  for (i = 0; i < m; i++)
    for(j = 0; j < n; j++)
    {
      A[n * i + j] = 0;
    }
}

//Matrix Multiplication Routine
// C = A*B
/*
void MatrixMultiply(float_prec* A, float_prec* B, int m, int p, int n, float_prec* C)
{
	// A = input matrix (m x p)
	// B = input matrix (p x n)
	// m = rows in A
	// p = columns in A = number of rows in B
	// n = columns in B
	// C = output matrix = A*B (m x n)
	int i, j, k;
	for (i = 0; i < m; i++)
		for(j = 0; j < n; j++)
		{
			C[n * i + j] = 0;
			for (k = 0; k < p; k++)
				C[n * i + j] = C[n * i + j] + A[p * i + k] * B[n * k + j];
		}
}*/

/*  Edited: scan by row first
 *  Rewrite: C <= A*B
 *
 *  col_A = row_B ?
 *
 *  row_C = row_A;
 *  col_C = col_B;
 */
int MatrixMultiply(float_prec* C, float_prec* A, int row_A, int col_A, float_prec* B, int row_B, int col_B){
	int index, row, col, row_C, col_C;
	if (col_A != row_B) { matrix_error_handle(); return 0;}
	row_C = row_A;
	col_C = col_B;
	for(row = 0; row < row_C; row++){
		for (col = 0; col < col_C; col++){
			C[col_C * row + col] = 0;
			for (index = 0; index < col_A; index++){
				C[col_C * row + col] = C[col_C * row + col] + A[col_A * row + index] * B[col_B * index + col];
			}
		}
	}
	return 1;
}



/*  Edited: scan by row first
 *  Rewrite: C = A*B
 *
 *  col_A = row_B ?
 *
 *  row_C = row_A;
 *  col_C = col_B;

float_prec* MatrixMultiply(float_prec* A, int row_A, int col_A, float_prec* B, int row_B, int col_B){
	int index, row, col;
	if (col_A != row_B) { matrix_error_handle(); return 0;}
	const int row_C = row_A;
	const int col_C = col_B;
	float_prec C[row_C*col_C];
	for(row = 0; row < row_C; row++){
		for (col = 0; col < col_C; col++){
			C[col_C * row + col] = 0;
			for (index = 0; index < col_A; index++){
				C[col_C * row + col] = C[col_C * row + col] + A[col_A * row + index] * B[col_B * index + col];
			}
		}
	}
	return C;
}
 */


/*	Edited: scan by row first
 *	C <= A & B
 *
 *  Axis 1: Horizontal
 *  row_A = row_B = row_C ?
 *  col_C = col_A + col_B;
 *
 *  Axis 0: Vertical
 *  row_C = row_A + row_B;
 *  col_A = col_B = col_C ?
 *
 *  // Axis 2: Join two array
 *  // row_A = row_B?
 *  // col_C = col_A + col_B;
 */
int MatrixConcatenate(float_prec* C, float_prec* A, int row_A, int col_A, float_prec* B, int row_B, int col_B, int axis){
	int col, row, col_C, row_C;
	if (axis > 2 || axis < 0) { matrix_error_handle(); return 0;}
	if (axis == 1){ // horizontal
		if (row_A != row_B) { matrix_error_handle(); return 0;}
		row_C = row_A;
		col_C = col_A + col_B;
		for(row = 0; row < row_C; row++){
			for (col = 0; col < col_C; col++){
				if (col > (col_A-1))
					C[col_C * row + col] = B[col_B * row + (col-col_A)];
				else
					C[col_C * row + col] = A[col_A * row + col];
			}
		}
	}
	if (axis == 0){ // vertical
		if (col_A != col_B) { matrix_error_handle(); return 0;}
		row_C = row_A + row_B;
		col_C = col_A;
		for(row = 0; row < row_C; row++){
			for (col = 0; col < col_C; col++){
				if (row > (row_A-1))
					C[col_C * row + col] = B[col_B * (row-row_A) + col];
				else
					C[col_C * row + col] = A[col_A * row + col];
			}
		}
	}
	/*
	if (axis == 2){ // Join
			if (row_A != row_B) { matrix_error_handle(); return 0;}
			if (row_A != 1) { matrix_error_handle(); return 0;}
			row_C = row_A;
			col_C = col_A + col_B;
			for(row = 0; row < row_C; row++){
				for (col = 0; col < col_C; col++){
					if (col > (col_A-1))
						C[col_C * row + col] = B[col_B * row + (col-col_A)];
					else
						C[col_C * row + col] = A[col_A * row + col];
				}
			}
		}
		*/
	return 1;
}

//Matrix Addition Routine
// C <= A + B
void MatrixAdd(float_prec* C, float_prec* A, float_prec* B, int m, int n){
	// A = input matrix (m x n)
	// B = input matrix (m x n)
	// m = number of rows in A = number of rows in B
	// n = number of columns in A = number of columns in B
	// C = output matrix = A+B (m x n)
	int i, j;
	for (i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			C[n * i + j] = A[n * i + j] + B[n * i + j];
}


//Matrix Subtraction Routine
// C <= A - B
void MatrixSubtract(float_prec* C, float_prec* A, float_prec* B, int m, int n){
	// A = input matrix (m x n)
	// B = input matrix (m x n)
	// m = rows in A = number of rows in B
	// n = columns in A = number of columns in B
	// C = output matrix = A-B (m x n)
	int i, j;
	for (i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			C[n * i + j] = A[n * i + j] - B[n * i + j];
}


//Matrix Transpose Routine
// C <= A
void MatrixTranspose(float_prec* C, float_prec* A, int m, int n){
	// A = input matrix (m x n)
	// m = number of rows in A
	// n = number of columns in A
	// C = output matrix = the transpose of A (n x m)
	int i, j;
	for (i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			C[m * j + i] = A[n * i + j];
}

// A <= A * k
void MatrixScale(float_prec* A, int m, int n, float_prec k){
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[n * i + j] = A[n * i + j] * k;
}


//Matrix Inversion Routine
// * This function inverts a matrix based on the Gauss Jordan method.
// * Specifically, it uses partial pivoting to improve numeric stability.
// * The algorithm is drawn from those presented in
//	 NUMERICAL RECIPES: The Art of Scientific Computing.
// * The function returns 1 on success, 0 on failure.
// * NOTE: The argument is ALSO the result matrix, meaning the input matrix is REPLACED
int MatrixInvert(float_prec* A, int n){
	// A = input matrix AND result matrix
	// n = number of rows = number of columns in A (n x n)
	int pivrow = 0;		// keeps track of current pivot row
	int k, i, j;		// k: overall index along diagonal; i: row index; j: col index
	int pivrows[n]; // keeps track of rows swaps to undo at end
	float_prec tmp;		// used for finding max value and making column swaps

	for (k = 0; k < n; k++)
	{
		// find pivot row, the row with biggest entry in current column
		tmp = 0;
		for (i = k; i < n; i++)
		{
			if (fabs(A[i * n + k]) >= tmp)	// 'Avoid using other functions inside abs()?'
			{
				tmp = fabs(A[i * n + k]);
				pivrow = i;
			}
		}

		// check for singular matrix
		if (A[pivrow * n + k] == 0.0f)
		{
			//Serial.println("Inversion failed due to singular matrix");
			matrix_error_handle_inv(1);
			return 0;
		}

		// Execute pivot (row swap) if needed
		if (pivrow != k)
		{
			// swap row k with pivrow
			for (j = 0; j < n; j++)
			{
				tmp = A[k * n + j];
				A[k * n + j] = A[pivrow * n + j];
				A[pivrow * n + j] = tmp;
			}
		}
		pivrows[k] = pivrow;	// record row swap (even if no swap happened)

		tmp = 1.0f / A[k * n + k];	// invert pivot element
		A[k * n + k] = 1.0f;		// This element of input matrix becomes result matrix

		// Perform row reduction (divide every element by pivot)
		for (j = 0; j < n; j++)
		{
			A[k * n + j] = A[k * n + j] * tmp;
		}

		// Now eliminate all other entries in this column
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				tmp = A[i * n + k];
				A[i * n + k] = 0.0f; // The other place where in matrix becomes result mat
				for (j = 0; j < n; j++)
				{
					A[i * n + j] = A[i * n + j] - A[k * n + j] * tmp;
				}
			}
		}
	}

	// Done, now need to undo pivot row swaps by doing column swaps in reverse order
	for (k = n - 1; k >= 0; k--)
	{
		if (pivrows[k] != k)
		{
			for (i = 0; i < n; i++)
			{
				tmp = A[i * n + k];
				A[i * n + k] = A[i * n + pivrows[k]];
				A[i * n + pivrows[k]] = tmp;
			}
		}
	}
	matrix_error_handle_inv(0);
	return 1;
}
