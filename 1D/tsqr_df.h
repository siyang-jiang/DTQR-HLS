#include "hls_stream.h"
#include "hls_math.h"
#include "ap_fixed.h"
#include "ap_int.h"

#define COLS 16
#define ROWS 2 * COLS
#define LEN (COLS + 1) * COLS / 2
#define Q_LEN (COLS + 1) * COLS
#include <iostream>

using namespace hls;

//typedef float MATRIX_T;

typedef ap_fixed<14,2,AP_RND_CONV> fixed_cs;
typedef ap_ufixed<26,1,AP_RND_CONV> fixed_cs_mag;
typedef ap_fixed<16,10,AP_RND_CONV> MATRIX_T;
typedef ap_ufixed<23,18> fixed_double;
typedef ap_ufixed<16,10> fixed_u;

typedef ap_uint<10> large_int;
typedef ap_uint<6> medium_int;
typedef ap_uint<2> small_int;


typedef int uint_i;
int top(MATRIX_T A1[LEN],
		MATRIX_T A2[LEN],
		fixed_cs Q_L[Q_LEN],
		fixed_cs Q_R[Q_LEN],
		MATRIX_T R[LEN]);

MATRIX_T qrf_mag(MATRIX_T a, MATRIX_T b);
void qrf_mm(fixed_cs c, fixed_cs s, MATRIX_T &op1, MATRIX_T &op2);

