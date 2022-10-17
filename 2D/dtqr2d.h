#ifndef qr_GR_H
#define qr_GR_H
#include "hls_stream.h"
#include "hls_math.h"
#include "ap_fixed.h"
#include "ap_int.h"
#define COLS 24
#define BUFFER 2
#define ROWS 2 * COLS
#define LEN (COLS + 1) * COLS / 2
#define Q_LEN (COLS + 1) * COLS
#define LEN_S (COLS) * (COLS - 1) / 2

#include <iostream>

using namespace hls;

typedef ap_fixed<14,2,AP_RND_CONV> fixed_cs;
typedef ap_ufixed<26,1,AP_RND_CONV> fixed_cs_mag;
typedef ap_fixed<16,10,AP_RND_CONV> MATRIX_T;
typedef ap_ufixed<23,18> fixed_double;
typedef ap_ufixed<16,10> fixed_u;

typedef ap_uint<10> large_int;
typedef ap_uint<6> medium_int;
typedef ap_uint<2> small_int;


typedef int uint_i;
int top(
	stream<MATRIX_T>&A1,
	stream<MATRIX_T>&A2,
	stream<fixed_cs>&Q_L,
	stream<fixed_cs>&Q_R,
	stream<MATRIX_T>&R
);

MATRIX_T qrf_mag(MATRIX_T a, MATRIX_T b);
void qrf_mm(MATRIX_T c, MATRIX_T s, MATRIX_T &op1, MATRIX_T &op2);

#endif
