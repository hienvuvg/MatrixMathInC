/*
 *  Created by Hien Vu on 04/21/2021.
 *  Based on MatrixMath.h for Arduino
 */

#ifndef _MATRIX_H_
#define _MATRIX_H_

#define HORIZONTAL 1
#define VERTICAL 0

/* Set this define to choose math precision of the system */
#define PRECISION_SINGLE    1
#define PRECISION_DOUBLE    2
#define FPU_PRECISION       (PRECISION_SINGLE)

#if (FPU_PRECISION == PRECISION_SINGLE)
    #define float_prec          float
    #define float_prec_ZERO     (1e-7)
    #define float_prec_ZERO_ECO (1e-5)      /* 'Economical' zero, for noisy calculation where 'somewhat zero' is good enough */
#elif (FPU_PRECISION == PRECISION_DOUBLE)
    #define float_prec          double
    #define float_prec_ZERO     (1e-13)
    #define float_prec_ZERO_ECO (1e-8)      /* 'Economical' zero, for noisy calculation where 'somewhat zero' is good enough */
#else
    #error("FPU_PRECISION has not been defined!");
#endif

	//MatrixMath();
	//void Print(float_prec* A, int m, int n, String label);
	void MatrixCopy(float_prec* A, float_prec* B, int row_A, int col_A);
	void MatrixWrite(float_prec* A, int row_A, int col_A, float_prec data[]);
	void MatrixIdentity(float_prec* A, int n);
	void MatrixZeros(float_prec* A, int n, int m);
	int MatrixConcatenate(float_prec* C, float_prec* A, int row_A, int col_A, float_prec* B, int row_B, int col_B, int axis);
	int MatrixMultiply(float_prec* C, float_prec* A, int row_A, int col_A, float_prec* B, int row_B, int col_B);
	//float_prec* Multiply(float_prec* A, int row_A, int col_A, float_prec* B, int row_B, int col_B);
	void MatrixAdd(float_prec* C, float_prec* A, float_prec* B, int m, int n);
	void MatrixSubtract(float_prec* C, float_prec* A, float_prec* B, int m, int n);
	void MatrixTranspose(float_prec* C, float_prec* A, int m, int n);
	void MatrixScale(float_prec* A, int m, int n, float_prec k);
	int MatrixInvert(float_prec* A, int n);

void matrix_error_handle (void);
void matrix_error_handle_inv (int led_on);

#endif
