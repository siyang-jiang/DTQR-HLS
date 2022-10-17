// TSQR 1D
#include "tsqr_df.h"
#include <iostream>

const int ROWS_C = ROWS;
const int COLS_C = COLS;
const int COLS_C_N = COLS_C - 1;
const int FIFO_LEN = LEN;
const fixed_cs_mag one = 1;

using namespace std;


MATRIX_T qrf_mag(MATRIX_T a, MATRIX_T b,fixed_cs &c, fixed_cs &s) // computes the magnitude of a and b
{
	fixed_double temp = a*a + b*b;;
	fixed_u mag = hls::sqrt(temp);
	fixed_cs_mag rmag = one/mag;
	c = a*rmag; 
	s = b*rmag;
	return mag;
}

void qrf_mm_Q(fixed_cs c, fixed_cs s, fixed_cs &op1, fixed_cs &op2)
{
#pragma HLS inline
	fixed_cs a = op2 * s + op1 * c; // multi
	fixed_cs b = op2 * c - op1 * s;
	op1 = a;
	op2 = b;
}

void qrf_mm(fixed_cs c, fixed_cs s, MATRIX_T &op1, MATRIX_T &op2)
{
#pragma HLS inline
	MATRIX_T a = op2 * s + op1 * c;
	MATRIX_T b = op2 * c - op1 * s;
	op1 = a;
	op2 = b;
}

void feeder(MATRIX_T A1[LEN],
			MATRIX_T A2[LEN],
			stream<MATRIX_T> &FeedDiagUp,
			stream<MATRIX_T> &FeedDiagLow,
			stream<MATRIX_T> &FeedElseUp,
			stream<MATRIX_T> &FeedElseLow)
{
	uint_i count = 0;
	for (uint_i j = 0; j < COLS; ++j)
	{
#pragma HLS loop_flatten off
		for (uint_i i = j; i < COLS; ++i)
		{
#pragma HLS pipeline II = 1
			if (i == j)
			{
				FeedDiagUp.write(A1[count]);
				FeedDiagLow.write(A2[count]);
				count += 1;
			}
			else
			{
				FeedElseUp.write(A1[count]);
				FeedElseLow.write(A2[count]);
				count += 1;
			}
		}
	}
}

void PE_head(stream<MATRIX_T> &in_A_diag_up,
			 stream<MATRIX_T> &in_A_diag_low,
			 stream<MATRIX_T> &in_A_else_up,
			 stream<MATRIX_T> &in_A_else_low,

			 stream<MATRIX_T> &out_A_diag_up,
			 stream<MATRIX_T> &out_A_else_up,
			 stream<MATRIX_T> &out_A_else_low,

			 stream<fixed_cs> &pass_c,
			 stream<fixed_cs> &pass_s,

			 stream<fixed_cs> &final_c,
			 stream<fixed_cs> &final_s,
			 stream<MATRIX_T> &out_R)
{
	fixed_cs c,s;
	MATRIX_T A_up[COLS];
#pragma HLS ARRAY_PARTITION variable = A_up complete dim = 1
	MATRIX_T A_low[COLS];
#pragma HLS ARRAY_PARTITION variable = A_low complete dim = 1
	for (uint_i i = 0; i < COLS; i++)
	{
#pragma HLS pipeline II = 1
		A_up[i] = in_A_diag_up.read();
		A_low[i] = in_A_diag_low.read();
	}

	for (uint_i j = 0; j < (LEN - COLS); j++)
	{
#pragma HLS pipeline II = 1
		out_A_else_up.write(in_A_else_up.read());
		out_A_else_low.write(in_A_else_low.read());
	}
	for (uint_i i = 0; i < COLS; i++)
	{
#pragma HLS unroll factor = COLS_C
		MATRIX_T mag = qrf_mag(A_up[i], A_low[i],c,s);
		A_up[i] = mag;
		A_low[i] = 0;
		final_c.write(c);
		final_s.write(s);
		pass_c.write(c);
		pass_s.write(s);
	}
	for (int i = 1; i < COLS; i++)
	{
#pragma HLS pipeline II = 1
		out_A_diag_up.write(A_up[i]);
	}
	out_R.write(A_up[0]);
}

void PE(uint_i id,

		stream<MATRIX_T> &in_A_diag_up,
		stream<MATRIX_T> &in_A_else_up,
		stream<MATRIX_T> &in_A_else_low,

		stream<MATRIX_T> &out_A_diag_up,
		stream<MATRIX_T> &out_A_else_up,
		stream<MATRIX_T> &out_A_else_low,

		stream<fixed_cs> &in_c,
		stream<fixed_cs> &in_s,

		stream<fixed_cs> &out_c,
		stream<fixed_cs> &out_s,

		stream<fixed_cs> &final_c,
		stream<fixed_cs> &final_s,

		stream<MATRIX_T> &out_R)
{

	fixed_cs C[COLS] = {0};
	fixed_cs S[COLS] = {0};
	MATRIX_T x, y;
Loop_1:
	for (int i = 0; i < COLS - id + 1; ++i)
	{
#pragma HLS pipeline II = 1
#pragma HLS loop_tripcount min = 0 max = COLS_C
		C[i] = in_c.read();
		S[i] = in_s.read();
	}
	uint_i count = (id - 1) * COLS - ((id - 1) * (id - 2)) / 2 + 1;
	stream<MATRIX_T> y_s;
#pragma HLS stream variable = y_s depth = FIFO_LEN
Loop_3:
	for (uint_i j = id; j < COLS; ++j)
	{
#pragma HLS loop_flatten off
#pragma HLS loop_tripcount min = 0 max = COLS_C
	Loop_3_1:
		for (uint_i i = j; i < COLS; ++i)
		{
#pragma HLS loop_tripcount min = 0 max = COLS_C
#pragma HLS pipeline II = 1
			x = in_A_else_up.read();
			y = in_A_else_low.read();

			qrf_mm(C[j - id], S[j - id], x, y);

			if (j == id)
			{
				out_R.write(x);
			}
			else
			{
				out_A_else_up.write(x);
			}
			y_s.write(y);
		}
	}
	
	MATRIX_T up[COLS], low[COLS];
Loop_4:
	for (uint_i i = 0; i < COLS - id; ++i)
	{
#pragma HLS pipeline II = 1
#pragma HLS loop_tripcount min = 1 max = COLS_C
		up[i] = in_A_diag_up.read();
	}
Loop_5:
	for (uint_i i = 0; i < COLS - id; ++i)
	{
#pragma HLS loop_tripcount min = 1 max = COLS_C
#pragma HLS loop_flatten off

Loop_6:
		for (uint_i j = i; j < COLS - id; ++j)
		{
#pragma HLS pipeline II = 1
#pragma HLS loop_tripcount min = 1 max = COLS_C
			if (i == j)
				low[i] = y_s.read();
			else
				out_A_else_low.write(y_s.read());
		}
	}
Loop_7:
	fixed_cs c,s;
	for (int i = 0; i < COLS - id; ++i)
	{
#pragma HLS unroll factor = COLS_C
#pragma HLS loop_tripcount min = 1 max = COLS_C
		MATRIX_T mag = qrf_mag(up[i], low[i],c,s);
		up[i] = mag;
		low[i] = 0;
		out_c.write(c);
		out_s.write(s);

		final_c.write(c);
		final_s.write(s);
	}
Loop_8:
	for (uint_i i = 1; i < COLS - id; i++)
	{
#pragma HLS pipeline II = 1
#pragma HLS loop_tripcount min = 2 max = COLS_C
		out_A_diag_up.write(up[i]);
	}
	out_R.write(up[0]);
}

void PE_tail(
	stream<MATRIX_T> &in_A_diag_up,
	stream<MATRIX_T> &in_A_else_up,
	stream<MATRIX_T> &in_A_else_low,

	stream<fixed_cs> &in_c,
	stream<fixed_cs> &in_s,

	stream<fixed_cs> &final_c,
	stream<fixed_cs> &final_s,

	stream<MATRIX_T> &out_R)
{

#pragma HLS pipeline II = 1

	MATRIX_T x, y, A_temp;
	fixed_cs c,s,c1,s1;
	c = in_c.read();
	s = in_s.read();
	x = in_A_else_up.read();
	y = in_A_else_low.read();
	A_temp = in_A_diag_up.read();

	qrf_mm(c, s, x, y);
	out_R.write(x);

	c = in_c.read();
	s = in_s.read();

	MATRIX_T mag = qrf_mag(A_temp, y,c1,s1);

	final_c.write(c1);
	final_s.write(s1);

	out_R.write(mag);
}

void PE3_head(stream<fixed_cs> &in_c,
			  stream<fixed_cs> &in_s,
			  stream<fixed_cs> &out_Q_left,
			  stream<fixed_cs> &out_Q_right,
			  stream<fixed_cs> &final_Q_left,
			  stream<fixed_cs> &final_Q_right)
{

	fixed_cs c, s;
	fixed_cs x, y;
	for (medium_int i = 0; i < COLS_C; i++)
	{
		c = in_c.read();
		s = in_s.read();
		for (small_int j = 0; j < 2; j++)
		{
#pragma HLS pipeline II = 1
			if (j==0){
				x = 1;
				y = 0;
			}else{
				x = 0;
				y = 1;
			}
			qrf_mm_Q(c, s, x, y);
			if (i == 0){
				final_Q_left.write(x);
				out_Q_right.write(y);
			}
			else if (i == COLS_C_N){
				out_Q_left.write(x);
				final_Q_right.write(y);
			}
			else{
				out_Q_left.write(x);
				out_Q_right.write(y);
			}
		}
	}
}

void PE3(medium_int id,
		 stream<fixed_cs> &in_c,
		 stream<fixed_cs> &in_s,
		 stream<fixed_cs> &in_Q_left,
		 stream<fixed_cs> &in_Q_right,
		 stream<fixed_cs> &out_Q_left,
		 stream<fixed_cs> &out_Q_right,
		 stream<fixed_cs> &final_Q_left,
		 stream<fixed_cs> &final_Q_right)
{
	fixed_cs c, s;
	fixed_cs x, y;
	for (medium_int i = 0; i < COLS; i++)
	{
		if (i < COLS-id){
		c = in_c.read();
		s = in_s.read();
		for (small_int _ = 0; _ < 2; _++)
		{
			for (medium_int j = 0; j < COLS; j++)
			{
				if (j < id+1)
				{
#pragma HLS pipeline II = 1
					if (j == 0)
						x = 0;
					else
						x = in_Q_left.read();

					if (j == id)
						y = 0;
					else
						y = in_Q_right.read();

					qrf_mm_Q(c, s, x, y);

					if (i == 0)
						final_Q_left.write(x);
					else
						out_Q_left.write(x);

					if (i == COLS - id - 1)
						final_Q_right.write(y);
					else
						out_Q_right.write(y);
				}
			}
		}
		}
	}
}

void PE3_tail(
	medium_int id, // id start at 0
	stream<fixed_cs> &in_c,
	stream<fixed_cs> &in_s,
	stream<fixed_cs> &in_Q_left,
	stream<fixed_cs> &in_Q_right,
	stream<fixed_cs> &final_Q_left,
	stream<fixed_cs> &final_Q_right)
{
	fixed_cs c, s;
	fixed_cs x, y;
	c = in_c.read();
	s = in_s.read();
	for (small_int _ = 0; _ < 2; _++)
	{
		x = 0;
		y = in_Q_right.read();
		qrf_mm_Q(c, s, x, y);
		final_Q_left.write(x);
		final_Q_right.write(y);
		for (medium_int j = 0; j < COLS-2; j++)
#pragma HLS pipeline II = 1
		{
			if (j < id-1)
			{
				x = in_Q_left.read();
				y = in_Q_right.read();
				qrf_mm_Q(c, s, x, y);
				final_Q_left.write(x);
				final_Q_right.write(y);
			}
		}
		x = in_Q_left.read();
		y = 0;
		qrf_mm_Q(c, s, x, y);
		final_Q_left.write(x);
		final_Q_right.write(y);
	}
}

void collector_Q(
	stream<fixed_cs> *final_Q_left,
	stream<fixed_cs> *final_Q_right,
	fixed_cs Q_L[Q_LEN],
	fixed_cs Q_R[Q_LEN])
{
	int count;
	for (int i = 0; i < COLS; ++i)
	{
#pragma HLS pipeline II = 1
		for (int j = 0; j < i + 1; ++j)
		{
#pragma HLS pipeline II = 1
			count = ((i * (i + 1)) / 2 + j) * 2;
			Q_L[count] = final_Q_left[i].read();
			Q_R[count] = final_Q_right[i].read();
			Q_L[count + 1] = final_Q_left[i].read();
			Q_R[count + 1] = final_Q_right[i].read();
		}
	}
}

void collector_R(stream<MATRIX_T> *in_R, MATRIX_T R[LEN])
{
	int counts = 0;
	R[0] = in_R[0].read();
	for (int i = 0; i < COLS - 1; i++)
	{
		for (int j = 0; j < COLS - i; j++)
		{
			counts = 1 + (2 * COLS - i + 1) * i / 2;
			R[counts + j] = in_R[i + 1].read();
		}
	}
}

int top(MATRIX_T A1[LEN],
		MATRIX_T A2[LEN],
		fixed_cs Q_L[Q_LEN],
		fixed_cs Q_R[Q_LEN],
		MATRIX_T R[LEN])
{
#pragma HLS dataflow

	stream<MATRIX_T> Feed[4];
	stream<fixed_cs> Q_left[COLS], Q_right[COLS];
	stream<fixed_cs> final_Q_left[COLS], final_Q_right[COLS];
	stream<MATRIX_T> out_A_diag_up[COLS_C_N];
	stream<MATRIX_T> out_A_else_up[COLS_C_N];
	stream<MATRIX_T> out_A_else_low[COLS_C_N];
	stream<fixed_cs> pass_c[COLS];
	stream<fixed_cs> pass_s[COLS];
	stream<MATRIX_T> out_R[COLS];
	stream<fixed_cs> final_c[COLS], final_s[COLS];

#pragma HLS resource variable = Q_left core = FIFO_BRAM
#pragma HLS resource variable = Q_right core = FIFO_BRAM
#pragma HLS resource variable = out_A_diag_up core = FIFO_BRAM
#pragma HLS resource variable = out_A_else_up core = FIFO_BRAM
#pragma HLS resource variable = out_A_else_low core = FIFO_BRAM
#pragma HLS resource variable = out_R core = FIFO_BRAM
#pragma HLS resource variable = pass_c core = FIFO_SRL
#pragma HLS resource variable = pass_s core = FIFO_SRL
#pragma HLS resource variable = final_c core = FIFO_BRAM
#pragma HLS resource variable = final_s core = FIFO_BRAM
#pragma HLS resource variable = final_Q_left core = FIFO_BRAM
#pragma HLS resource variable = final_Q_right core = FIFO_BRAM

	feeder(A1, A2,
		   Feed[0],
		   Feed[1],
		   Feed[2],
		   Feed[3]);

	PE_head(Feed[0], Feed[1], Feed[2], Feed[3],
			out_A_diag_up[0], out_A_else_up[0], out_A_else_low[0],
			pass_c[0], pass_s[0],
			final_c[0], final_s[0],
			out_R[0]);

	PE(1, out_A_diag_up[0], out_A_else_up[0], out_A_else_low[0],
	   out_A_diag_up[1], out_A_else_up[1], out_A_else_low[1],
	   pass_c[0], pass_s[0], pass_c[1], pass_s[1],
	   final_c[1], final_s[1],
	   out_R[1]);

	PE(2, out_A_diag_up[1], out_A_else_up[1], out_A_else_low[1],
	   out_A_diag_up[2], out_A_else_up[2], out_A_else_low[2],
	   pass_c[1], pass_s[1], pass_c[2], pass_s[2],
	   final_c[2], final_s[2],
	   out_R[2]);

	PE(3, out_A_diag_up[2], out_A_else_up[2], out_A_else_low[2],
	   out_A_diag_up[3], out_A_else_up[3], out_A_else_low[3],
	   pass_c[2], pass_s[2], pass_c[3], pass_s[3],
	   final_c[3], final_s[3],
	   out_R[3]);

	PE(4, out_A_diag_up[3], out_A_else_up[3], out_A_else_low[3],
	   out_A_diag_up[4], out_A_else_up[4], out_A_else_low[4],
	   pass_c[3], pass_s[3], pass_c[4], pass_s[4],
	   final_c[4], final_s[4],
	   out_R[4]);

	PE(5, out_A_diag_up[4], out_A_else_up[4], out_A_else_low[4],
	   out_A_diag_up[5], out_A_else_up[5], out_A_else_low[5],
	   pass_c[4], pass_s[4], pass_c[5], pass_s[5],
	   final_c[5], final_s[5],
	   out_R[5]);

	PE(6, out_A_diag_up[5], out_A_else_up[5], out_A_else_low[5],
	   out_A_diag_up[6], out_A_else_up[6], out_A_else_low[6],
	   pass_c[5], pass_s[5], pass_c[6], pass_s[6],
	   final_c[6], final_s[6],
	   out_R[6]);

	PE(7, out_A_diag_up[6], out_A_else_up[6], out_A_else_low[6],
	   out_A_diag_up[7], out_A_else_up[7], out_A_else_low[7],
	   pass_c[6], pass_s[6], pass_c[7], pass_s[7],
	   final_c[7], final_s[7],
	   out_R[7]);

	PE(8, out_A_diag_up[7], out_A_else_up[7], out_A_else_low[7],
	   out_A_diag_up[8], out_A_else_up[8], out_A_else_low[8],
	   pass_c[7], pass_s[7], pass_c[8], pass_s[8],
	   final_c[8], final_s[8],
	   out_R[8]);

	PE(9, out_A_diag_up[8], out_A_else_up[8], out_A_else_low[8],
	   out_A_diag_up[9], out_A_else_up[9], out_A_else_low[9],
	   pass_c[8], pass_s[8], pass_c[9], pass_s[9],
	   final_c[9], final_s[9],
	   out_R[9]);

	PE(10, out_A_diag_up[9], out_A_else_up[9], out_A_else_low[9],
	   out_A_diag_up[10], out_A_else_up[10], out_A_else_low[10],
	   pass_c[9], pass_s[9], pass_c[10], pass_s[10],
	   final_c[10], final_s[10],
	   out_R[10]);

	PE(11, out_A_diag_up[10], out_A_else_up[10], out_A_else_low[10],
	   out_A_diag_up[11], out_A_else_up[11], out_A_else_low[11],
	   pass_c[10], pass_s[10], pass_c[11], pass_s[11],
	   final_c[11], final_s[11],
	   out_R[11]);

	PE(12, out_A_diag_up[11], out_A_else_up[11], out_A_else_low[11],
	   out_A_diag_up[12], out_A_else_up[12], out_A_else_low[12],
	   pass_c[11], pass_s[11], pass_c[12], pass_s[12],
	   final_c[12], final_s[12],
	   out_R[12]);

	PE(13, out_A_diag_up[12], out_A_else_up[12], out_A_else_low[12],
	   out_A_diag_up[13], out_A_else_up[13], out_A_else_low[13],
	   pass_c[12], pass_s[12], pass_c[13], pass_s[13],
	   final_c[13], final_s[13],
	   out_R[13]);

	PE(14, out_A_diag_up[13], out_A_else_up[13], out_A_else_low[13],
	   out_A_diag_up[14], out_A_else_up[14], out_A_else_low[14],
	   pass_c[13], pass_s[13], pass_c[14], pass_s[14],
	   final_c[14], final_s[14],
	   out_R[14]);

//	PE(15, out_A_diag_up[14], out_A_else_up[14], out_A_else_low[14],
//	   out_A_diag_up[15], out_A_else_up[15], out_A_else_low[15],
//	   pass_c[14], pass_s[14], pass_c[15], pass_s[15],
//	   final_c[15], final_s[15],
//	   out_R[15]);

//	PE(16, out_A_diag_up[15], out_A_else_up[15], out_A_else_low[15],
//	   out_A_diag_up[16], out_A_else_up[16], out_A_else_low[16],
//	   pass_c[15], pass_s[15], pass_c[16], pass_s[16],
//	   final_c[16], final_s[16],
//	   out_R[16]);
//
//	PE(17, out_A_diag_up[16], out_A_else_up[16], out_A_else_low[16],
//	   out_A_diag_up[17], out_A_else_up[17], out_A_else_low[17],
//	   pass_c[16], pass_s[16], pass_c[17], pass_s[17],
//	   final_c[17], final_s[17],
//	   out_R[17]);
//
//	PE(18, out_A_diag_up[17], out_A_else_up[17], out_A_else_low[17],
//	   out_A_diag_up[18], out_A_else_up[18], out_A_else_low[18],
//	   pass_c[17], pass_s[17], pass_c[18], pass_s[18],
//	   final_c[18], final_s[18],
//	   out_R[18]);
//
//	PE(19, out_A_diag_up[18], out_A_else_up[18], out_A_else_low[18],
//	   out_A_diag_up[19], out_A_else_up[19], out_A_else_low[19],
//	   pass_c[18], pass_s[18], pass_c[19], pass_s[19],
//	   final_c[19], final_s[19],
//	   out_R[19]);
//
//	PE(20, out_A_diag_up[19], out_A_else_up[19], out_A_else_low[19],
//	   out_A_diag_up[20], out_A_else_up[20], out_A_else_low[20],
//	   pass_c[19], pass_s[19], pass_c[20], pass_s[20],
//	   final_c[20], final_s[20],
//	   out_R[20]);
//
//	PE(21, out_A_diag_up[20], out_A_else_up[20], out_A_else_low[20],
//	   out_A_diag_up[21], out_A_else_up[21], out_A_else_low[21],
//	   pass_c[20], pass_s[20], pass_c[21], pass_s[21],
//	   final_c[21], final_s[21],
//	   out_R[21]);
//
//	PE(22, out_A_diag_up[21], out_A_else_up[21], out_A_else_low[21],
//	   out_A_diag_up[22], out_A_else_up[22], out_A_else_low[22],
//	   pass_c[21], pass_s[21], pass_c[22], pass_s[22],
//	   final_c[22], final_s[22],
//	   out_R[22]);

	PE_tail(out_A_diag_up[COLS - 2], out_A_else_up[COLS - 2], out_A_else_low[COLS - 2],
			pass_c[COLS - 2], pass_s[COLS - 2], final_c[COLS - 1], final_s[COLS - 1], out_R[COLS - 1]);

	collector_R(out_R, R);

	PE3_head(final_c[0], final_s[0], Q_left[0], Q_right[0], final_Q_left[0], final_Q_right[0]);
	PE3(1, final_c[1], final_s[1], Q_left[0], Q_right[0], Q_left[1], Q_right[1], final_Q_left[1], final_Q_right[1]);
	PE3(2, final_c[2], final_s[2], Q_left[1], Q_right[1], Q_left[2], Q_right[2], final_Q_left[2], final_Q_right[2]);
	PE3(3, final_c[3], final_s[3], Q_left[2], Q_right[2], Q_left[3], Q_right[3], final_Q_left[3], final_Q_right[3]);
	PE3(4, final_c[4], final_s[4], Q_left[3], Q_right[3], Q_left[4], Q_right[4], final_Q_left[4], final_Q_right[4]);
	PE3(5, final_c[5], final_s[5], Q_left[4], Q_right[4], Q_left[5], Q_right[5], final_Q_left[5], final_Q_right[5]);
	PE3(6, final_c[6], final_s[6], Q_left[5], Q_right[5], Q_left[6], Q_right[6], final_Q_left[6], final_Q_right[6]);
	PE3(7, final_c[7], final_s[7], Q_left[6], Q_right[6], Q_left[7], Q_right[7], final_Q_left[7], final_Q_right[7]);
	PE3(8, final_c[8], final_s[8], Q_left[7], Q_right[7], Q_left[8], Q_right[8], final_Q_left[8], final_Q_right[8]);
	PE3(9, final_c[9], final_s[9], Q_left[8], Q_right[8], Q_left[9], Q_right[9], final_Q_left[9], final_Q_right[9]);
	PE3(10, final_c[10], final_s[10], Q_left[9], Q_right[9], Q_left[10], Q_right[10], final_Q_left[10], final_Q_right[10]);
	PE3(11, final_c[11], final_s[11], Q_left[10], Q_right[10], Q_left[11], Q_right[11], final_Q_left[11], final_Q_right[11]);
	PE3(12, final_c[12], final_s[12], Q_left[11], Q_right[11], Q_left[12], Q_right[12], final_Q_left[12], final_Q_right[12]);
	PE3(13, final_c[13], final_s[13], Q_left[12], Q_right[12], Q_left[13], Q_right[13], final_Q_left[13], final_Q_right[13]);
	PE3(14, final_c[14], final_s[14], Q_left[13], Q_right[13], Q_left[14], Q_right[14], final_Q_left[14], final_Q_right[14]);
//	PE3(15, final_c[15], final_s[15], Q_left[14], Q_right[14], Q_left[15], Q_right[15], final_Q_left[15], final_Q_right[15]);
//	PE3(16, final_c[16], final_s[16], Q_left[15], Q_right[15], Q_left[16], Q_right[16], final_Q_left[16], final_Q_right[16]);
//	PE3(17, final_c[17], final_s[17], Q_left[16], Q_right[16], Q_left[17], Q_right[17], final_Q_left[17], final_Q_right[17]);
//	PE3(18, final_c[18], final_s[18], Q_left[17], Q_right[17], Q_left[18], Q_right[18], final_Q_left[18], final_Q_right[18]);
//	PE3(19, final_c[19], final_s[19], Q_left[18], Q_right[18], Q_left[19], Q_right[19], final_Q_left[19], final_Q_right[19]);
//	PE3(20, final_c[20], final_s[20], Q_left[19], Q_right[19], Q_left[20], Q_right[20], final_Q_left[20], final_Q_right[20]);
//	PE3(21, final_c[21], final_s[21], Q_left[20], Q_right[20], Q_left[21], Q_right[21], final_Q_left[21], final_Q_right[21]);
//	PE3(22, final_c[22], final_s[22], Q_left[21], Q_right[21], Q_left[22], Q_right[22], final_Q_left[22], final_Q_right[22]);
	PE3_tail(COLS - 1, final_c[COLS - 1], final_s[COLS - 1], Q_left[COLS - 2], Q_right[COLS - 2], final_Q_left[COLS - 1], final_Q_right[COLS - 1]);

	collector_Q(final_Q_left, final_Q_right, Q_L, Q_R);

	return 0;
}
