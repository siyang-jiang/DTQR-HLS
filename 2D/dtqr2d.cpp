#include "dtqr2d.h"

#include <iostream>

const int ROWS_C = ROWS;
const int COLS_C = COLS;
const int COLS_C_N = COLS_C - 1;
const int LEN_S_C = LEN_S;
const fixed_cs_mag one = 1;
using namespace std;

MATRIX_T qrf_mag_B(MATRIX_T a, MATRIX_T b,fixed_cs &c, fixed_cs &s)
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

// Feeder to prepare the data;

void feeder(
	stream<MATRIX_T> &A1,
	stream<MATRIX_T> &A2,
	stream<MATRIX_T> *Feed_A,
	stream<MATRIX_T> *Feed_B)
{
	medium_int count = 0;
	medium_int j = 0;
	medium_int max = COLS;
	for (large_int sum = 0; sum < LEN; sum++)
	{
#pragma HLS pipeline II = 1
		Feed_A[count].write(A1.read());
		Feed_B[count].write(A2.read());
		count++;
		if (count == max){
			count = 0;
			max --;
		}
	}
}

void PE1(
	medium_int output_size, // id start at 0
	stream<MATRIX_T> &in_a,
	stream<MATRIX_T> &in_b,
	stream<MATRIX_T> &out_a,
	stream<fixed_cs> &out_c,
	stream<fixed_cs> &out_s,
	stream<fixed_cs> &out_c_Q,
	stream<fixed_cs> &out_s_Q,
	stream<MATRIX_T> &out_R)
{
#pragma HLS FUNCTION_INSTANTIATE variable=output_size
	MATRIX_T A, B;
	fixed_cs c[COLS], s[COLS]; 
	MATRIX_T mag;
	int p = 0;
loading_Loop:
	for (int i = 0; i < COLS_C; i++)
	{
#pragma HLS loop_tripcount min=1 max=COLS_C avg=COLS_C/2
#pragma HLS pipeline II = 3
		A = in_a.read();
		B = in_b.read();
		mag = qrf_mag_B(A, B,c[i],s[i]);
		if (i == 0)
			out_R.write(mag);
		else
			out_a.write(mag);
		if (i == output_size-1)
			break;
	}

transfer_Loop_1:
	for (int i = 0; i < COLS_C; i++) 
	{
#pragma HLS pipeline II = 1
		out_c_Q.write(c[i]);
		out_s_Q.write(s[i]);
		if (i != output_size-1){
			out_c.write(c[i]);
			out_s.write(s[i]);
		}else{
			break;
		}
	}
}

void PE1_tail(
	stream<MATRIX_T> &in_a,
	stream<MATRIX_T> &in_b,
	stream<fixed_cs> &out_c_Q,
	stream<fixed_cs> &out_s_Q,
	stream<MATRIX_T> &out_R)
{

	MATRIX_T A = in_a.read();
	MATRIX_T B = in_b.read();
	fixed_cs c,s;
	MATRIX_T mag = qrf_mag_B(A, B,c,s);
	out_c_Q.write(c);
	out_s_Q.write(s);
	out_R.write(mag);
}

void PE2(
	medium_int output_size, //output size = COLS - idx - idy - 1
	stream<MATRIX_T> &in_a,
	stream<MATRIX_T> &in_b,

	stream<fixed_cs> &in_c,
	stream<fixed_cs> &in_s,

	stream<MATRIX_T> &pass_a,
	stream<MATRIX_T> &pass_b,

	stream<fixed_cs> &pass_c,
	stream<fixed_cs> &pass_s,

	stream<MATRIX_T> &final_R)

{
	MATRIX_T A, B;
	fixed_cs c, s;
	medium_int pass_cs_range = output_size-1;
	for (medium_int i = 0; i < COLS_C; i++)
	{
#pragma HLS pipeline II = 5
		if (i < output_size){
			c = in_c.read();
			s = in_s.read();
			if (i != pass_cs_range){
				pass_c.write(c);
				pass_s.write(s);
			}
			A = in_a.read();
			B= in_b.read();
			qrf_mm(c, s, A, B);
			if (i == 0)
				final_R.write(A);
			else
				pass_a.write(A);
			pass_b.write(B);
		}
	}
}

void PE2_tail(
	stream<MATRIX_T> &in_a,
	stream<MATRIX_T> &in_b,
	stream<fixed_cs> &in_c,
	stream<fixed_cs> &in_s,
	stream<MATRIX_T> &out_b,
	stream<MATRIX_T> &out_a_R)
{

	MATRIX_T A, B;
	fixed_cs c, s;

	c = in_c.read();
	s = in_s.read();
	A = in_a.read();
	B = in_b.read();
	qrf_mm(c, s, A, B);
	out_a_R.write(A);
	out_b.write(B);
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
// #pragma HLS FUNCTION_INSTANTIATE variable=id
	fixed_cs c, s;
	fixed_cs x, y;
	for (medium_int i = 0; i < COLS; i++)
	{
		if (i < COLS-id){
		c = in_c.read();
		s = in_s.read();
	// }
	// for (medium_int i = 0; i < COLS - id; i++){
		for (small_int _ = 0; _ < 2; _++)
		{
			for (medium_int j = 0; j < COLS; j++)
			{
#pragma HLS pipeline II = 1
				if (j < id+1)
				{
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
	stream<fixed_cs> &Q_L,
	stream<fixed_cs> &Q_R)
{
	large_int count = 0;
	for (medium_int i = 0; i < COLS; ++i)
	{
		for (medium_int j = 0; j < COLS; ++j)
		{
			if (j < i+1)
			{
				for (small_int _ = 0; _<2; _++){
					#pragma HLS pipeline II = 1
					Q_L.write(final_Q_left[i].read());
					Q_R.write(final_Q_right[i].read());
					count++;
				}
			}
		}
	}
}

void collector_R(stream<MATRIX_T> *in_R, stream<MATRIX_T> &R)
{
	for (int i = 0; i < LEN; i++)
	{
		R.write(in_R[i].read());
	}
}

void region(
	stream<MATRIX_T,COLS> *Feed_A,
	stream<MATRIX_T,COLS> *Feed_B,
	stream<fixed_cs,COLS*2> *final_Q_left,
	stream<fixed_cs,COLS*2> *final_Q_right,
	stream<MATRIX_T,2> *final_R
){
	#pragma HLS dataflow
	#pragma HLS INTERFACE mode=ap_ctrl_none port=return

	stream<MATRIX_T,COLS> out_a[COLS]; 
	stream<fixed_cs,COLS> PE1_pass_c[COLS], PE1_pass_s[COLS];
	stream<fixed_cs,2> pass_c[COLS][COLS], pass_s[COLS][COLS];  
	stream<MATRIX_T,COLS> pass_a[LEN_S], pass_b[LEN_S];
	stream<fixed_cs,COLS> out_c_Q[COLS], out_s_Q[COLS];
	stream<fixed_cs,COLS*2> Q_left[COLS], Q_right[COLS]; 

// REPLACE START
	PE1(24, Feed_A[0], Feed_B[0], out_a[0], PE1_pass_c[0], PE1_pass_s[0], out_c_Q[0], out_s_Q[0], final_R[0]);

	PE3_head(out_c_Q[0], out_s_Q[0], Q_left[0], Q_right[0], final_Q_left[0], final_Q_right[0]);
	PE2(23, Feed_A[1], Feed_B[1], PE1_pass_c[0], PE1_pass_s[0], pass_a[0], pass_b[0], pass_c[0][0], pass_s[0][0], final_R[1]);

	PE1(23, out_a[0], pass_b[0], out_a[1], PE1_pass_c[1], PE1_pass_s[1], out_c_Q[1], out_s_Q[1], final_R[24]);
	PE2(22, Feed_A[2], Feed_B[2], pass_c[0][0], pass_s[0][0], pass_a[1], pass_b[1], pass_c[0][1], pass_s[0][1], final_R[2]);

	PE3(1,out_c_Q[1], out_s_Q[1],Q_left[0], Q_right[0],Q_left[1], Q_right[1], final_Q_left[1], final_Q_right[1]);
	PE2(22, pass_a[0], pass_b[1], PE1_pass_c[1], PE1_pass_s[1], pass_a[22], pass_b[23], pass_c[1][0], pass_s[1][0], final_R[25]);
	PE2(21, Feed_A[3], Feed_B[3], pass_c[0][1], pass_s[0][1], pass_a[2], pass_b[2], pass_c[0][2], pass_s[0][2], final_R[3]);

	PE1(22, out_a[1], pass_b[23], out_a[2], PE1_pass_c[2], PE1_pass_s[2], out_c_Q[2], out_s_Q[2], final_R[47]);
	PE2(21, pass_a[1], pass_b[2], pass_c[1][0], pass_s[1][0], pass_a[23], pass_b[24], pass_c[1][1], pass_s[1][1], final_R[26]);
	PE2(20, Feed_A[4], Feed_B[4], pass_c[0][2], pass_s[0][2], pass_a[3], pass_b[3], pass_c[0][3], pass_s[0][3], final_R[4]);

	PE3(2,out_c_Q[2], out_s_Q[2],Q_left[1], Q_right[1],Q_left[2], Q_right[2], final_Q_left[2], final_Q_right[2]);
	PE2(21, pass_a[22], pass_b[24], PE1_pass_c[2], PE1_pass_s[2], pass_a[43], pass_b[45], pass_c[2][0], pass_s[2][0], final_R[48]);
	PE2(20, pass_a[2], pass_b[3], pass_c[1][1], pass_s[1][1], pass_a[24], pass_b[25], pass_c[1][2], pass_s[1][2], final_R[27]);
	PE2(19, Feed_A[5], Feed_B[5], pass_c[0][3], pass_s[0][3], pass_a[4], pass_b[4], pass_c[0][4], pass_s[0][4], final_R[5]);

	PE1(21, out_a[2], pass_b[45], out_a[3], PE1_pass_c[3], PE1_pass_s[3], out_c_Q[3], out_s_Q[3], final_R[69]);
	PE2(20, pass_a[23], pass_b[25], pass_c[2][0], pass_s[2][0], pass_a[44], pass_b[46], pass_c[2][1], pass_s[2][1], final_R[49]);
	PE2(19, pass_a[3], pass_b[4], pass_c[1][2], pass_s[1][2], pass_a[25], pass_b[26], pass_c[1][3], pass_s[1][3], final_R[28]);
	PE2(18, Feed_A[6], Feed_B[6], pass_c[0][4], pass_s[0][4], pass_a[5], pass_b[5], pass_c[0][5], pass_s[0][5], final_R[6]);

	PE3(3,out_c_Q[3], out_s_Q[3],Q_left[2], Q_right[2],Q_left[3], Q_right[3], final_Q_left[3], final_Q_right[3]);
	PE2(20, pass_a[43], pass_b[46], PE1_pass_c[3], PE1_pass_s[3], pass_a[63], pass_b[66], pass_c[3][0], pass_s[3][0], final_R[70]);
	PE2(19, pass_a[24], pass_b[26], pass_c[2][1], pass_s[2][1], pass_a[45], pass_b[47], pass_c[2][2], pass_s[2][2], final_R[50]);
	PE2(18, pass_a[4], pass_b[5], pass_c[1][3], pass_s[1][3], pass_a[26], pass_b[27], pass_c[1][4], pass_s[1][4], final_R[29]);
	PE2(17, Feed_A[7], Feed_B[7], pass_c[0][5], pass_s[0][5], pass_a[6], pass_b[6], pass_c[0][6], pass_s[0][6], final_R[7]);

	PE1(20, out_a[3], pass_b[66], out_a[4], PE1_pass_c[4], PE1_pass_s[4], out_c_Q[4], out_s_Q[4], final_R[90]);
	PE2(19, pass_a[44], pass_b[47], pass_c[3][0], pass_s[3][0], pass_a[64], pass_b[67], pass_c[3][1], pass_s[3][1], final_R[71]);
	PE2(18, pass_a[25], pass_b[27], pass_c[2][2], pass_s[2][2], pass_a[46], pass_b[48], pass_c[2][3], pass_s[2][3], final_R[51]);
	PE2(17, pass_a[5], pass_b[6], pass_c[1][4], pass_s[1][4], pass_a[27], pass_b[28], pass_c[1][5], pass_s[1][5], final_R[30]);
	PE2(16, Feed_A[8], Feed_B[8], pass_c[0][6], pass_s[0][6], pass_a[7], pass_b[7], pass_c[0][7], pass_s[0][7], final_R[8]);

	PE3(4,out_c_Q[4], out_s_Q[4],Q_left[3], Q_right[3],Q_left[4], Q_right[4], final_Q_left[4], final_Q_right[4]);
	PE2(19, pass_a[63], pass_b[67], PE1_pass_c[4], PE1_pass_s[4], pass_a[82], pass_b[86], pass_c[4][0], pass_s[4][0], final_R[91]);
	PE2(18, pass_a[45], pass_b[48], pass_c[3][1], pass_s[3][1], pass_a[65], pass_b[68], pass_c[3][2], pass_s[3][2], final_R[72]);
	PE2(17, pass_a[26], pass_b[28], pass_c[2][3], pass_s[2][3], pass_a[47], pass_b[49], pass_c[2][4], pass_s[2][4], final_R[52]);
	PE2(16, pass_a[6], pass_b[7], pass_c[1][5], pass_s[1][5], pass_a[28], pass_b[29], pass_c[1][6], pass_s[1][6], final_R[31]);
	PE2(15, Feed_A[9], Feed_B[9], pass_c[0][7], pass_s[0][7], pass_a[8], pass_b[8], pass_c[0][8], pass_s[0][8], final_R[9]);

	PE1(19, out_a[4], pass_b[86], out_a[5], PE1_pass_c[5], PE1_pass_s[5], out_c_Q[5], out_s_Q[5], final_R[110]);
	PE2(18, pass_a[64], pass_b[68], pass_c[4][0], pass_s[4][0], pass_a[83], pass_b[87], pass_c[4][1], pass_s[4][1], final_R[92]);
	PE2(17, pass_a[46], pass_b[49], pass_c[3][2], pass_s[3][2], pass_a[66], pass_b[69], pass_c[3][3], pass_s[3][3], final_R[73]);
	PE2(16, pass_a[27], pass_b[29], pass_c[2][4], pass_s[2][4], pass_a[48], pass_b[50], pass_c[2][5], pass_s[2][5], final_R[53]);
	PE2(15, pass_a[7], pass_b[8], pass_c[1][6], pass_s[1][6], pass_a[29], pass_b[30], pass_c[1][7], pass_s[1][7], final_R[32]);
	PE2(14, Feed_A[10], Feed_B[10], pass_c[0][8], pass_s[0][8], pass_a[9], pass_b[9], pass_c[0][9], pass_s[0][9], final_R[10]);

	PE3(5,out_c_Q[5], out_s_Q[5],Q_left[4], Q_right[4],Q_left[5], Q_right[5], final_Q_left[5], final_Q_right[5]);
	PE2(18, pass_a[82], pass_b[87], PE1_pass_c[5], PE1_pass_s[5], pass_a[100], pass_b[105], pass_c[5][0], pass_s[5][0], final_R[111]);
	PE2(17, pass_a[65], pass_b[69], pass_c[4][1], pass_s[4][1], pass_a[84], pass_b[88], pass_c[4][2], pass_s[4][2], final_R[93]);
	PE2(16, pass_a[47], pass_b[50], pass_c[3][3], pass_s[3][3], pass_a[67], pass_b[70], pass_c[3][4], pass_s[3][4], final_R[74]);
	PE2(15, pass_a[28], pass_b[30], pass_c[2][5], pass_s[2][5], pass_a[49], pass_b[51], pass_c[2][6], pass_s[2][6], final_R[54]);
	PE2(14, pass_a[8], pass_b[9], pass_c[1][7], pass_s[1][7], pass_a[30], pass_b[31], pass_c[1][8], pass_s[1][8], final_R[33]);
	PE2(13, Feed_A[11], Feed_B[11], pass_c[0][9], pass_s[0][9], pass_a[10], pass_b[10], pass_c[0][10], pass_s[0][10], final_R[11]);

	PE1(18, out_a[5], pass_b[105], out_a[6], PE1_pass_c[6], PE1_pass_s[6], out_c_Q[6], out_s_Q[6], final_R[129]);
	PE2(17, pass_a[83], pass_b[88], pass_c[5][0], pass_s[5][0], pass_a[101], pass_b[106], pass_c[5][1], pass_s[5][1], final_R[112]);
	PE2(16, pass_a[66], pass_b[70], pass_c[4][2], pass_s[4][2], pass_a[85], pass_b[89], pass_c[4][3], pass_s[4][3], final_R[94]);
	PE2(15, pass_a[48], pass_b[51], pass_c[3][4], pass_s[3][4], pass_a[68], pass_b[71], pass_c[3][5], pass_s[3][5], final_R[75]);
	PE2(14, pass_a[29], pass_b[31], pass_c[2][6], pass_s[2][6], pass_a[50], pass_b[52], pass_c[2][7], pass_s[2][7], final_R[55]);
	PE2(13, pass_a[9], pass_b[10], pass_c[1][8], pass_s[1][8], pass_a[31], pass_b[32], pass_c[1][9], pass_s[1][9], final_R[34]);
	PE2(12, Feed_A[12], Feed_B[12], pass_c[0][10], pass_s[0][10], pass_a[11], pass_b[11], pass_c[0][11], pass_s[0][11], final_R[12]);

	PE3(6,out_c_Q[6], out_s_Q[6],Q_left[5], Q_right[5],Q_left[6], Q_right[6], final_Q_left[6], final_Q_right[6]);
	PE2(17, pass_a[100], pass_b[106], PE1_pass_c[6], PE1_pass_s[6], pass_a[117], pass_b[123], pass_c[6][0], pass_s[6][0], final_R[130]);
	PE2(16, pass_a[84], pass_b[89], pass_c[5][1], pass_s[5][1], pass_a[102], pass_b[107], pass_c[5][2], pass_s[5][2], final_R[113]);
	PE2(15, pass_a[67], pass_b[71], pass_c[4][3], pass_s[4][3], pass_a[86], pass_b[90], pass_c[4][4], pass_s[4][4], final_R[95]);
	PE2(14, pass_a[49], pass_b[52], pass_c[3][5], pass_s[3][5], pass_a[69], pass_b[72], pass_c[3][6], pass_s[3][6], final_R[76]);
	PE2(13, pass_a[30], pass_b[32], pass_c[2][7], pass_s[2][7], pass_a[51], pass_b[53], pass_c[2][8], pass_s[2][8], final_R[56]);
	PE2(12, pass_a[10], pass_b[11], pass_c[1][9], pass_s[1][9], pass_a[32], pass_b[33], pass_c[1][10], pass_s[1][10], final_R[35]);
	PE2(11, Feed_A[13], Feed_B[13], pass_c[0][11], pass_s[0][11], pass_a[12], pass_b[12], pass_c[0][12], pass_s[0][12], final_R[13]);

	PE1(17, out_a[6], pass_b[123], out_a[7], PE1_pass_c[7], PE1_pass_s[7], out_c_Q[7], out_s_Q[7], final_R[147]);
	PE2(16, pass_a[101], pass_b[107], pass_c[6][0], pass_s[6][0], pass_a[118], pass_b[124], pass_c[6][1], pass_s[6][1], final_R[131]);
	PE2(15, pass_a[85], pass_b[90], pass_c[5][2], pass_s[5][2], pass_a[103], pass_b[108], pass_c[5][3], pass_s[5][3], final_R[114]);
	PE2(14, pass_a[68], pass_b[72], pass_c[4][4], pass_s[4][4], pass_a[87], pass_b[91], pass_c[4][5], pass_s[4][5], final_R[96]);
	PE2(13, pass_a[50], pass_b[53], pass_c[3][6], pass_s[3][6], pass_a[70], pass_b[73], pass_c[3][7], pass_s[3][7], final_R[77]);
	PE2(12, pass_a[31], pass_b[33], pass_c[2][8], pass_s[2][8], pass_a[52], pass_b[54], pass_c[2][9], pass_s[2][9], final_R[57]);
	PE2(11, pass_a[11], pass_b[12], pass_c[1][10], pass_s[1][10], pass_a[33], pass_b[34], pass_c[1][11], pass_s[1][11], final_R[36]);
	PE2(10, Feed_A[14], Feed_B[14], pass_c[0][12], pass_s[0][12], pass_a[13], pass_b[13], pass_c[0][13], pass_s[0][13], final_R[14]);

	PE3(7,out_c_Q[7], out_s_Q[7],Q_left[6], Q_right[6],Q_left[7], Q_right[7], final_Q_left[7], final_Q_right[7]);
	PE2(16, pass_a[117], pass_b[124], PE1_pass_c[7], PE1_pass_s[7], pass_a[133], pass_b[140], pass_c[7][0], pass_s[7][0], final_R[148]);
	PE2(15, pass_a[102], pass_b[108], pass_c[6][1], pass_s[6][1], pass_a[119], pass_b[125], pass_c[6][2], pass_s[6][2], final_R[132]);
	PE2(14, pass_a[86], pass_b[91], pass_c[5][3], pass_s[5][3], pass_a[104], pass_b[109], pass_c[5][4], pass_s[5][4], final_R[115]);
	PE2(13, pass_a[69], pass_b[73], pass_c[4][5], pass_s[4][5], pass_a[88], pass_b[92], pass_c[4][6], pass_s[4][6], final_R[97]);
	PE2(12, pass_a[51], pass_b[54], pass_c[3][7], pass_s[3][7], pass_a[71], pass_b[74], pass_c[3][8], pass_s[3][8], final_R[78]);
	PE2(11, pass_a[32], pass_b[34], pass_c[2][9], pass_s[2][9], pass_a[53], pass_b[55], pass_c[2][10], pass_s[2][10], final_R[58]);
	PE2(10, pass_a[12], pass_b[13], pass_c[1][11], pass_s[1][11], pass_a[34], pass_b[35], pass_c[1][12], pass_s[1][12], final_R[37]);
	PE2(9, Feed_A[15], Feed_B[15], pass_c[0][13], pass_s[0][13], pass_a[14], pass_b[14], pass_c[0][14], pass_s[0][14], final_R[15]);

	PE1(16, out_a[7], pass_b[140], out_a[8], PE1_pass_c[8], PE1_pass_s[8], out_c_Q[8], out_s_Q[8], final_R[164]);
	PE2(15, pass_a[118], pass_b[125], pass_c[7][0], pass_s[7][0], pass_a[134], pass_b[141], pass_c[7][1], pass_s[7][1], final_R[149]);
	PE2(14, pass_a[103], pass_b[109], pass_c[6][2], pass_s[6][2], pass_a[120], pass_b[126], pass_c[6][3], pass_s[6][3], final_R[133]);
	PE2(13, pass_a[87], pass_b[92], pass_c[5][4], pass_s[5][4], pass_a[105], pass_b[110], pass_c[5][5], pass_s[5][5], final_R[116]);
	PE2(12, pass_a[70], pass_b[74], pass_c[4][6], pass_s[4][6], pass_a[89], pass_b[93], pass_c[4][7], pass_s[4][7], final_R[98]);
	PE2(11, pass_a[52], pass_b[55], pass_c[3][8], pass_s[3][8], pass_a[72], pass_b[75], pass_c[3][9], pass_s[3][9], final_R[79]);
	PE2(10, pass_a[33], pass_b[35], pass_c[2][10], pass_s[2][10], pass_a[54], pass_b[56], pass_c[2][11], pass_s[2][11], final_R[59]);
	PE2(9, pass_a[13], pass_b[14], pass_c[1][12], pass_s[1][12], pass_a[35], pass_b[36], pass_c[1][13], pass_s[1][13], final_R[38]);
	PE2(8, Feed_A[16], Feed_B[16], pass_c[0][14], pass_s[0][14], pass_a[15], pass_b[15], pass_c[0][15], pass_s[0][15], final_R[16]);

	PE3(8,out_c_Q[8], out_s_Q[8],Q_left[7], Q_right[7],Q_left[8], Q_right[8], final_Q_left[8], final_Q_right[8]);
	PE2(15, pass_a[133], pass_b[141], PE1_pass_c[8], PE1_pass_s[8], pass_a[148], pass_b[156], pass_c[8][0], pass_s[8][0], final_R[165]);
	PE2(14, pass_a[119], pass_b[126], pass_c[7][1], pass_s[7][1], pass_a[135], pass_b[142], pass_c[7][2], pass_s[7][2], final_R[150]);
	PE2(13, pass_a[104], pass_b[110], pass_c[6][3], pass_s[6][3], pass_a[121], pass_b[127], pass_c[6][4], pass_s[6][4], final_R[134]);
	PE2(12, pass_a[88], pass_b[93], pass_c[5][5], pass_s[5][5], pass_a[106], pass_b[111], pass_c[5][6], pass_s[5][6], final_R[117]);
	PE2(11, pass_a[71], pass_b[75], pass_c[4][7], pass_s[4][7], pass_a[90], pass_b[94], pass_c[4][8], pass_s[4][8], final_R[99]);
	PE2(10, pass_a[53], pass_b[56], pass_c[3][9], pass_s[3][9], pass_a[73], pass_b[76], pass_c[3][10], pass_s[3][10], final_R[80]);
	PE2(9, pass_a[34], pass_b[36], pass_c[2][11], pass_s[2][11], pass_a[55], pass_b[57], pass_c[2][12], pass_s[2][12], final_R[60]);
	PE2(8, pass_a[14], pass_b[15], pass_c[1][13], pass_s[1][13], pass_a[36], pass_b[37], pass_c[1][14], pass_s[1][14], final_R[39]);
	PE2(7, Feed_A[17], Feed_B[17], pass_c[0][15], pass_s[0][15], pass_a[16], pass_b[16], pass_c[0][16], pass_s[0][16], final_R[17]);

	PE1(15, out_a[8], pass_b[156], out_a[9], PE1_pass_c[9], PE1_pass_s[9], out_c_Q[9], out_s_Q[9], final_R[180]);
	PE2(14, pass_a[134], pass_b[142], pass_c[8][0], pass_s[8][0], pass_a[149], pass_b[157], pass_c[8][1], pass_s[8][1], final_R[166]);
	PE2(13, pass_a[120], pass_b[127], pass_c[7][2], pass_s[7][2], pass_a[136], pass_b[143], pass_c[7][3], pass_s[7][3], final_R[151]);
	PE2(12, pass_a[105], pass_b[111], pass_c[6][4], pass_s[6][4], pass_a[122], pass_b[128], pass_c[6][5], pass_s[6][5], final_R[135]);
	PE2(11, pass_a[89], pass_b[94], pass_c[5][6], pass_s[5][6], pass_a[107], pass_b[112], pass_c[5][7], pass_s[5][7], final_R[118]);
	PE2(10, pass_a[72], pass_b[76], pass_c[4][8], pass_s[4][8], pass_a[91], pass_b[95], pass_c[4][9], pass_s[4][9], final_R[100]);
	PE2(9, pass_a[54], pass_b[57], pass_c[3][10], pass_s[3][10], pass_a[74], pass_b[77], pass_c[3][11], pass_s[3][11], final_R[81]);
	PE2(8, pass_a[35], pass_b[37], pass_c[2][12], pass_s[2][12], pass_a[56], pass_b[58], pass_c[2][13], pass_s[2][13], final_R[61]);
	PE2(7, pass_a[15], pass_b[16], pass_c[1][14], pass_s[1][14], pass_a[37], pass_b[38], pass_c[1][15], pass_s[1][15], final_R[40]);
	PE2(6, Feed_A[18], Feed_B[18], pass_c[0][16], pass_s[0][16], pass_a[17], pass_b[17], pass_c[0][17], pass_s[0][17], final_R[18]);

	PE3(9,out_c_Q[9], out_s_Q[9],Q_left[8], Q_right[8],Q_left[9], Q_right[9], final_Q_left[9], final_Q_right[9]);
	PE2(14, pass_a[148], pass_b[157], PE1_pass_c[9], PE1_pass_s[9], pass_a[162], pass_b[171], pass_c[9][0], pass_s[9][0], final_R[181]);
	PE2(13, pass_a[135], pass_b[143], pass_c[8][1], pass_s[8][1], pass_a[150], pass_b[158], pass_c[8][2], pass_s[8][2], final_R[167]);
	PE2(12, pass_a[121], pass_b[128], pass_c[7][3], pass_s[7][3], pass_a[137], pass_b[144], pass_c[7][4], pass_s[7][4], final_R[152]);
	PE2(11, pass_a[106], pass_b[112], pass_c[6][5], pass_s[6][5], pass_a[123], pass_b[129], pass_c[6][6], pass_s[6][6], final_R[136]);
	PE2(10, pass_a[90], pass_b[95], pass_c[5][7], pass_s[5][7], pass_a[108], pass_b[113], pass_c[5][8], pass_s[5][8], final_R[119]);
	PE2(9, pass_a[73], pass_b[77], pass_c[4][9], pass_s[4][9], pass_a[92], pass_b[96], pass_c[4][10], pass_s[4][10], final_R[101]);
	PE2(8, pass_a[55], pass_b[58], pass_c[3][11], pass_s[3][11], pass_a[75], pass_b[78], pass_c[3][12], pass_s[3][12], final_R[82]);
	PE2(7, pass_a[36], pass_b[38], pass_c[2][13], pass_s[2][13], pass_a[57], pass_b[59], pass_c[2][14], pass_s[2][14], final_R[62]);
	PE2(6, pass_a[16], pass_b[17], pass_c[1][15], pass_s[1][15], pass_a[38], pass_b[39], pass_c[1][16], pass_s[1][16], final_R[41]);
	PE2(5, Feed_A[19], Feed_B[19], pass_c[0][17], pass_s[0][17], pass_a[18], pass_b[18], pass_c[0][18], pass_s[0][18], final_R[19]);

	PE1(14, out_a[9], pass_b[171], out_a[10], PE1_pass_c[10], PE1_pass_s[10], out_c_Q[10], out_s_Q[10], final_R[195]);
	PE2(13, pass_a[149], pass_b[158], pass_c[9][0], pass_s[9][0], pass_a[163], pass_b[172], pass_c[9][1], pass_s[9][1], final_R[182]);
	PE2(12, pass_a[136], pass_b[144], pass_c[8][2], pass_s[8][2], pass_a[151], pass_b[159], pass_c[8][3], pass_s[8][3], final_R[168]);
	PE2(11, pass_a[122], pass_b[129], pass_c[7][4], pass_s[7][4], pass_a[138], pass_b[145], pass_c[7][5], pass_s[7][5], final_R[153]);
	PE2(10, pass_a[107], pass_b[113], pass_c[6][6], pass_s[6][6], pass_a[124], pass_b[130], pass_c[6][7], pass_s[6][7], final_R[137]);
	PE2(9, pass_a[91], pass_b[96], pass_c[5][8], pass_s[5][8], pass_a[109], pass_b[114], pass_c[5][9], pass_s[5][9], final_R[120]);
	PE2(8, pass_a[74], pass_b[78], pass_c[4][10], pass_s[4][10], pass_a[93], pass_b[97], pass_c[4][11], pass_s[4][11], final_R[102]);
	PE2(7, pass_a[56], pass_b[59], pass_c[3][12], pass_s[3][12], pass_a[76], pass_b[79], pass_c[3][13], pass_s[3][13], final_R[83]);
	PE2(6, pass_a[37], pass_b[39], pass_c[2][14], pass_s[2][14], pass_a[58], pass_b[60], pass_c[2][15], pass_s[2][15], final_R[63]);
	PE2(5, pass_a[17], pass_b[18], pass_c[1][16], pass_s[1][16], pass_a[39], pass_b[40], pass_c[1][17], pass_s[1][17], final_R[42]);
	PE2(4, Feed_A[20], Feed_B[20], pass_c[0][18], pass_s[0][18], pass_a[19], pass_b[19], pass_c[0][19], pass_s[0][19], final_R[20]);

	PE3(10,out_c_Q[10], out_s_Q[10],Q_left[9], Q_right[9],Q_left[10], Q_right[10], final_Q_left[10], final_Q_right[10]);
	PE2(13, pass_a[162], pass_b[172], PE1_pass_c[10], PE1_pass_s[10], pass_a[175], pass_b[185], pass_c[10][0], pass_s[10][0], final_R[196]);
	PE2(12, pass_a[150], pass_b[159], pass_c[9][1], pass_s[9][1], pass_a[164], pass_b[173], pass_c[9][2], pass_s[9][2], final_R[183]);
	PE2(11, pass_a[137], pass_b[145], pass_c[8][3], pass_s[8][3], pass_a[152], pass_b[160], pass_c[8][4], pass_s[8][4], final_R[169]);
	PE2(10, pass_a[123], pass_b[130], pass_c[7][5], pass_s[7][5], pass_a[139], pass_b[146], pass_c[7][6], pass_s[7][6], final_R[154]);
	PE2(9, pass_a[108], pass_b[114], pass_c[6][7], pass_s[6][7], pass_a[125], pass_b[131], pass_c[6][8], pass_s[6][8], final_R[138]);
	PE2(8, pass_a[92], pass_b[97], pass_c[5][9], pass_s[5][9], pass_a[110], pass_b[115], pass_c[5][10], pass_s[5][10], final_R[121]);
	PE2(7, pass_a[75], pass_b[79], pass_c[4][11], pass_s[4][11], pass_a[94], pass_b[98], pass_c[4][12], pass_s[4][12], final_R[103]);
	PE2(6, pass_a[57], pass_b[60], pass_c[3][13], pass_s[3][13], pass_a[77], pass_b[80], pass_c[3][14], pass_s[3][14], final_R[84]);
	PE2(5, pass_a[38], pass_b[40], pass_c[2][15], pass_s[2][15], pass_a[59], pass_b[61], pass_c[2][16], pass_s[2][16], final_R[64]);
	PE2(4, pass_a[18], pass_b[19], pass_c[1][17], pass_s[1][17], pass_a[40], pass_b[41], pass_c[1][18], pass_s[1][18], final_R[43]);
	PE2(3, Feed_A[21], Feed_B[21], pass_c[0][19], pass_s[0][19], pass_a[20], pass_b[20], pass_c[0][20], pass_s[0][20], final_R[21]);

	PE1(13, out_a[10], pass_b[185], out_a[11], PE1_pass_c[11], PE1_pass_s[11], out_c_Q[11], out_s_Q[11], final_R[209]);
	PE2(12, pass_a[163], pass_b[173], pass_c[10][0], pass_s[10][0], pass_a[176], pass_b[186], pass_c[10][1], pass_s[10][1], final_R[197]);
	PE2(11, pass_a[151], pass_b[160], pass_c[9][2], pass_s[9][2], pass_a[165], pass_b[174], pass_c[9][3], pass_s[9][3], final_R[184]);
	PE2(10, pass_a[138], pass_b[146], pass_c[8][4], pass_s[8][4], pass_a[153], pass_b[161], pass_c[8][5], pass_s[8][5], final_R[170]);
	PE2(9, pass_a[124], pass_b[131], pass_c[7][6], pass_s[7][6], pass_a[140], pass_b[147], pass_c[7][7], pass_s[7][7], final_R[155]);
	PE2(8, pass_a[109], pass_b[115], pass_c[6][8], pass_s[6][8], pass_a[126], pass_b[132], pass_c[6][9], pass_s[6][9], final_R[139]);
	PE2(7, pass_a[93], pass_b[98], pass_c[5][10], pass_s[5][10], pass_a[111], pass_b[116], pass_c[5][11], pass_s[5][11], final_R[122]);
	PE2(6, pass_a[76], pass_b[80], pass_c[4][12], pass_s[4][12], pass_a[95], pass_b[99], pass_c[4][13], pass_s[4][13], final_R[104]);
	PE2(5, pass_a[58], pass_b[61], pass_c[3][14], pass_s[3][14], pass_a[78], pass_b[81], pass_c[3][15], pass_s[3][15], final_R[85]);
	PE2(4, pass_a[39], pass_b[41], pass_c[2][16], pass_s[2][16], pass_a[60], pass_b[62], pass_c[2][17], pass_s[2][17], final_R[65]);
	PE2(3, pass_a[19], pass_b[20], pass_c[1][18], pass_s[1][18], pass_a[41], pass_b[42], pass_c[1][19], pass_s[1][19], final_R[44]);
	PE2(2, Feed_A[22], Feed_B[22], pass_c[0][20], pass_s[0][20], pass_a[21], pass_b[21], pass_c[0][21], pass_s[0][21], final_R[22]);

	PE3(11,out_c_Q[11], out_s_Q[11],Q_left[10], Q_right[10],Q_left[11], Q_right[11], final_Q_left[11], final_Q_right[11]);
	PE2(12, pass_a[175], pass_b[186], PE1_pass_c[11], PE1_pass_s[11], pass_a[187], pass_b[198], pass_c[11][0], pass_s[11][0], final_R[210]);
	PE2(11, pass_a[164], pass_b[174], pass_c[10][1], pass_s[10][1], pass_a[177], pass_b[187], pass_c[10][2], pass_s[10][2], final_R[198]);
	PE2(10, pass_a[152], pass_b[161], pass_c[9][3], pass_s[9][3], pass_a[166], pass_b[175], pass_c[9][4], pass_s[9][4], final_R[185]);
	PE2(9, pass_a[139], pass_b[147], pass_c[8][5], pass_s[8][5], pass_a[154], pass_b[162], pass_c[8][6], pass_s[8][6], final_R[171]);
	PE2(8, pass_a[125], pass_b[132], pass_c[7][7], pass_s[7][7], pass_a[141], pass_b[148], pass_c[7][8], pass_s[7][8], final_R[156]);
	PE2(7, pass_a[110], pass_b[116], pass_c[6][9], pass_s[6][9], pass_a[127], pass_b[133], pass_c[6][10], pass_s[6][10], final_R[140]);
	PE2(6, pass_a[94], pass_b[99], pass_c[5][11], pass_s[5][11], pass_a[112], pass_b[117], pass_c[5][12], pass_s[5][12], final_R[123]);
	PE2(5, pass_a[77], pass_b[81], pass_c[4][13], pass_s[4][13], pass_a[96], pass_b[100], pass_c[4][14], pass_s[4][14], final_R[105]);
	PE2(4, pass_a[59], pass_b[62], pass_c[3][15], pass_s[3][15], pass_a[79], pass_b[82], pass_c[3][16], pass_s[3][16], final_R[86]);
	PE2(3, pass_a[40], pass_b[42], pass_c[2][17], pass_s[2][17], pass_a[61], pass_b[63], pass_c[2][18], pass_s[2][18], final_R[66]);
	PE2(2, pass_a[20], pass_b[21], pass_c[1][19], pass_s[1][19], pass_a[42], pass_b[43], pass_c[1][20], pass_s[1][20], final_R[45]);
	PE2_tail(Feed_A[23], Feed_B[23], pass_c[0][21], pass_s[0][21], pass_b[22], final_R[23]);

	PE1(12, out_a[11], pass_b[198], out_a[12], PE1_pass_c[12], PE1_pass_s[12], out_c_Q[12], out_s_Q[12], final_R[222]);
	PE2(11, pass_a[176], pass_b[187], pass_c[11][0], pass_s[11][0], pass_a[188], pass_b[199], pass_c[11][1], pass_s[11][1], final_R[211]);
	PE2(10, pass_a[165], pass_b[175], pass_c[10][2], pass_s[10][2], pass_a[178], pass_b[188], pass_c[10][3], pass_s[10][3], final_R[199]);
	PE2(9, pass_a[153], pass_b[162], pass_c[9][4], pass_s[9][4], pass_a[167], pass_b[176], pass_c[9][5], pass_s[9][5], final_R[186]);
	PE2(8, pass_a[140], pass_b[148], pass_c[8][6], pass_s[8][6], pass_a[155], pass_b[163], pass_c[8][7], pass_s[8][7], final_R[172]);
	PE2(7, pass_a[126], pass_b[133], pass_c[7][8], pass_s[7][8], pass_a[142], pass_b[149], pass_c[7][9], pass_s[7][9], final_R[157]);
	PE2(6, pass_a[111], pass_b[117], pass_c[6][10], pass_s[6][10], pass_a[128], pass_b[134], pass_c[6][11], pass_s[6][11], final_R[141]);
	PE2(5, pass_a[95], pass_b[100], pass_c[5][12], pass_s[5][12], pass_a[113], pass_b[118], pass_c[5][13], pass_s[5][13], final_R[124]);
	PE2(4, pass_a[78], pass_b[82], pass_c[4][14], pass_s[4][14], pass_a[97], pass_b[101], pass_c[4][15], pass_s[4][15], final_R[106]);
	PE2(3, pass_a[60], pass_b[63], pass_c[3][16], pass_s[3][16], pass_a[80], pass_b[83], pass_c[3][17], pass_s[3][17], final_R[87]);
	PE2(2, pass_a[41], pass_b[43], pass_c[2][18], pass_s[2][18], pass_a[62], pass_b[64], pass_c[2][19], pass_s[2][19], final_R[67]);
	PE2_tail(pass_a[21], pass_b[22], pass_c[1][20], pass_s[1][20], pass_b[44], final_R[46]);

	PE3(12,out_c_Q[12], out_s_Q[12],Q_left[11], Q_right[11],Q_left[12], Q_right[12], final_Q_left[12], final_Q_right[12]);
	PE2(11, pass_a[187], pass_b[199], PE1_pass_c[12], PE1_pass_s[12], pass_a[198], pass_b[210], pass_c[12][0], pass_s[12][0], final_R[223]);
	PE2(10, pass_a[177], pass_b[188], pass_c[11][1], pass_s[11][1], pass_a[189], pass_b[200], pass_c[11][2], pass_s[11][2], final_R[212]);
	PE2(9, pass_a[166], pass_b[176], pass_c[10][3], pass_s[10][3], pass_a[179], pass_b[189], pass_c[10][4], pass_s[10][4], final_R[200]);
	PE2(8, pass_a[154], pass_b[163], pass_c[9][5], pass_s[9][5], pass_a[168], pass_b[177], pass_c[9][6], pass_s[9][6], final_R[187]);
	PE2(7, pass_a[141], pass_b[149], pass_c[8][7], pass_s[8][7], pass_a[156], pass_b[164], pass_c[8][8], pass_s[8][8], final_R[173]);
	PE2(6, pass_a[127], pass_b[134], pass_c[7][9], pass_s[7][9], pass_a[143], pass_b[150], pass_c[7][10], pass_s[7][10], final_R[158]);
	PE2(5, pass_a[112], pass_b[118], pass_c[6][11], pass_s[6][11], pass_a[129], pass_b[135], pass_c[6][12], pass_s[6][12], final_R[142]);
	PE2(4, pass_a[96], pass_b[101], pass_c[5][13], pass_s[5][13], pass_a[114], pass_b[119], pass_c[5][14], pass_s[5][14], final_R[125]);
	PE2(3, pass_a[79], pass_b[83], pass_c[4][15], pass_s[4][15], pass_a[98], pass_b[102], pass_c[4][16], pass_s[4][16], final_R[107]);
	PE2(2, pass_a[61], pass_b[64], pass_c[3][17], pass_s[3][17], pass_a[81], pass_b[84], pass_c[3][18], pass_s[3][18], final_R[88]);
	PE2_tail(pass_a[42], pass_b[44], pass_c[2][19], pass_s[2][19], pass_b[65], final_R[68]);

	PE1(11, out_a[12], pass_b[210], out_a[13], PE1_pass_c[13], PE1_pass_s[13], out_c_Q[13], out_s_Q[13], final_R[234]);
	PE2(10, pass_a[188], pass_b[200], pass_c[12][0], pass_s[12][0], pass_a[199], pass_b[211], pass_c[12][1], pass_s[12][1], final_R[224]);
	PE2(9, pass_a[178], pass_b[189], pass_c[11][2], pass_s[11][2], pass_a[190], pass_b[201], pass_c[11][3], pass_s[11][3], final_R[213]);
	PE2(8, pass_a[167], pass_b[177], pass_c[10][4], pass_s[10][4], pass_a[180], pass_b[190], pass_c[10][5], pass_s[10][5], final_R[201]);
	PE2(7, pass_a[155], pass_b[164], pass_c[9][6], pass_s[9][6], pass_a[169], pass_b[178], pass_c[9][7], pass_s[9][7], final_R[188]);
	PE2(6, pass_a[142], pass_b[150], pass_c[8][8], pass_s[8][8], pass_a[157], pass_b[165], pass_c[8][9], pass_s[8][9], final_R[174]);
	PE2(5, pass_a[128], pass_b[135], pass_c[7][10], pass_s[7][10], pass_a[144], pass_b[151], pass_c[7][11], pass_s[7][11], final_R[159]);
	PE2(4, pass_a[113], pass_b[119], pass_c[6][12], pass_s[6][12], pass_a[130], pass_b[136], pass_c[6][13], pass_s[6][13], final_R[143]);
	PE2(3, pass_a[97], pass_b[102], pass_c[5][14], pass_s[5][14], pass_a[115], pass_b[120], pass_c[5][15], pass_s[5][15], final_R[126]);
	PE2(2, pass_a[80], pass_b[84], pass_c[4][16], pass_s[4][16], pass_a[99], pass_b[103], pass_c[4][17], pass_s[4][17], final_R[108]);
	PE2_tail(pass_a[62], pass_b[65], pass_c[3][18], pass_s[3][18], pass_b[85], final_R[89]);

	PE3(13,out_c_Q[13], out_s_Q[13],Q_left[12], Q_right[12],Q_left[13], Q_right[13], final_Q_left[13], final_Q_right[13]);
	PE2(10, pass_a[198], pass_b[211], PE1_pass_c[13], PE1_pass_s[13], pass_a[208], pass_b[221], pass_c[13][0], pass_s[13][0], final_R[235]);
	PE2(9, pass_a[189], pass_b[201], pass_c[12][1], pass_s[12][1], pass_a[200], pass_b[212], pass_c[12][2], pass_s[12][2], final_R[225]);
	PE2(8, pass_a[179], pass_b[190], pass_c[11][3], pass_s[11][3], pass_a[191], pass_b[202], pass_c[11][4], pass_s[11][4], final_R[214]);
	PE2(7, pass_a[168], pass_b[178], pass_c[10][5], pass_s[10][5], pass_a[181], pass_b[191], pass_c[10][6], pass_s[10][6], final_R[202]);
	PE2(6, pass_a[156], pass_b[165], pass_c[9][7], pass_s[9][7], pass_a[170], pass_b[179], pass_c[9][8], pass_s[9][8], final_R[189]);
	PE2(5, pass_a[143], pass_b[151], pass_c[8][9], pass_s[8][9], pass_a[158], pass_b[166], pass_c[8][10], pass_s[8][10], final_R[175]);
	PE2(4, pass_a[129], pass_b[136], pass_c[7][11], pass_s[7][11], pass_a[145], pass_b[152], pass_c[7][12], pass_s[7][12], final_R[160]);
	PE2(3, pass_a[114], pass_b[120], pass_c[6][13], pass_s[6][13], pass_a[131], pass_b[137], pass_c[6][14], pass_s[6][14], final_R[144]);
	PE2(2, pass_a[98], pass_b[103], pass_c[5][15], pass_s[5][15], pass_a[116], pass_b[121], pass_c[5][16], pass_s[5][16], final_R[127]);
	PE2_tail(pass_a[81], pass_b[85], pass_c[4][17], pass_s[4][17], pass_b[104], final_R[109]);

	PE1(10, out_a[13], pass_b[221], out_a[14], PE1_pass_c[14], PE1_pass_s[14], out_c_Q[14], out_s_Q[14], final_R[245]);
	PE2(9, pass_a[199], pass_b[212], pass_c[13][0], pass_s[13][0], pass_a[209], pass_b[222], pass_c[13][1], pass_s[13][1], final_R[236]);
	PE2(8, pass_a[190], pass_b[202], pass_c[12][2], pass_s[12][2], pass_a[201], pass_b[213], pass_c[12][3], pass_s[12][3], final_R[226]);
	PE2(7, pass_a[180], pass_b[191], pass_c[11][4], pass_s[11][4], pass_a[192], pass_b[203], pass_c[11][5], pass_s[11][5], final_R[215]);
	PE2(6, pass_a[169], pass_b[179], pass_c[10][6], pass_s[10][6], pass_a[182], pass_b[192], pass_c[10][7], pass_s[10][7], final_R[203]);
	PE2(5, pass_a[157], pass_b[166], pass_c[9][8], pass_s[9][8], pass_a[171], pass_b[180], pass_c[9][9], pass_s[9][9], final_R[190]);
	PE2(4, pass_a[144], pass_b[152], pass_c[8][10], pass_s[8][10], pass_a[159], pass_b[167], pass_c[8][11], pass_s[8][11], final_R[176]);
	PE2(3, pass_a[130], pass_b[137], pass_c[7][12], pass_s[7][12], pass_a[146], pass_b[153], pass_c[7][13], pass_s[7][13], final_R[161]);
	PE2(2, pass_a[115], pass_b[121], pass_c[6][14], pass_s[6][14], pass_a[132], pass_b[138], pass_c[6][15], pass_s[6][15], final_R[145]);
	PE2_tail(pass_a[99], pass_b[104], pass_c[5][16], pass_s[5][16], pass_b[122], final_R[128]);

	PE3(14,out_c_Q[14], out_s_Q[14],Q_left[13], Q_right[13],Q_left[14], Q_right[14], final_Q_left[14], final_Q_right[14]);
	PE2(9, pass_a[208], pass_b[222], PE1_pass_c[14], PE1_pass_s[14], pass_a[217], pass_b[231], pass_c[14][0], pass_s[14][0], final_R[246]);
	PE2(8, pass_a[200], pass_b[213], pass_c[13][1], pass_s[13][1], pass_a[210], pass_b[223], pass_c[13][2], pass_s[13][2], final_R[237]);
	PE2(7, pass_a[191], pass_b[203], pass_c[12][3], pass_s[12][3], pass_a[202], pass_b[214], pass_c[12][4], pass_s[12][4], final_R[227]);
	PE2(6, pass_a[181], pass_b[192], pass_c[11][5], pass_s[11][5], pass_a[193], pass_b[204], pass_c[11][6], pass_s[11][6], final_R[216]);
	PE2(5, pass_a[170], pass_b[180], pass_c[10][7], pass_s[10][7], pass_a[183], pass_b[193], pass_c[10][8], pass_s[10][8], final_R[204]);
	PE2(4, pass_a[158], pass_b[167], pass_c[9][9], pass_s[9][9], pass_a[172], pass_b[181], pass_c[9][10], pass_s[9][10], final_R[191]);
	PE2(3, pass_a[145], pass_b[153], pass_c[8][11], pass_s[8][11], pass_a[160], pass_b[168], pass_c[8][12], pass_s[8][12], final_R[177]);
	PE2(2, pass_a[131], pass_b[138], pass_c[7][13], pass_s[7][13], pass_a[147], pass_b[154], pass_c[7][14], pass_s[7][14], final_R[162]);
	PE2_tail(pass_a[116], pass_b[122], pass_c[6][15], pass_s[6][15], pass_b[139], final_R[146]);

	PE1(9, out_a[14], pass_b[231], out_a[15], PE1_pass_c[15], PE1_pass_s[15], out_c_Q[15], out_s_Q[15], final_R[255]);
	PE2(8, pass_a[209], pass_b[223], pass_c[14][0], pass_s[14][0], pass_a[218], pass_b[232], pass_c[14][1], pass_s[14][1], final_R[247]);
	PE2(7, pass_a[201], pass_b[214], pass_c[13][2], pass_s[13][2], pass_a[211], pass_b[224], pass_c[13][3], pass_s[13][3], final_R[238]);
	PE2(6, pass_a[192], pass_b[204], pass_c[12][4], pass_s[12][4], pass_a[203], pass_b[215], pass_c[12][5], pass_s[12][5], final_R[228]);
	PE2(5, pass_a[182], pass_b[193], pass_c[11][6], pass_s[11][6], pass_a[194], pass_b[205], pass_c[11][7], pass_s[11][7], final_R[217]);
	PE2(4, pass_a[171], pass_b[181], pass_c[10][8], pass_s[10][8], pass_a[184], pass_b[194], pass_c[10][9], pass_s[10][9], final_R[205]);
	PE2(3, pass_a[159], pass_b[168], pass_c[9][10], pass_s[9][10], pass_a[173], pass_b[182], pass_c[9][11], pass_s[9][11], final_R[192]);
	PE2(2, pass_a[146], pass_b[154], pass_c[8][12], pass_s[8][12], pass_a[161], pass_b[169], pass_c[8][13], pass_s[8][13], final_R[178]);
	PE2_tail(pass_a[132], pass_b[139], pass_c[7][14], pass_s[7][14], pass_b[155], final_R[163]);

	PE3(15,out_c_Q[15], out_s_Q[15],Q_left[14], Q_right[14],Q_left[15], Q_right[15], final_Q_left[15], final_Q_right[15]);
	PE2(8, pass_a[217], pass_b[232], PE1_pass_c[15], PE1_pass_s[15], pass_a[225], pass_b[240], pass_c[15][0], pass_s[15][0], final_R[256]);
	PE2(7, pass_a[210], pass_b[224], pass_c[14][1], pass_s[14][1], pass_a[219], pass_b[233], pass_c[14][2], pass_s[14][2], final_R[248]);
	PE2(6, pass_a[202], pass_b[215], pass_c[13][3], pass_s[13][3], pass_a[212], pass_b[225], pass_c[13][4], pass_s[13][4], final_R[239]);
	PE2(5, pass_a[193], pass_b[205], pass_c[12][5], pass_s[12][5], pass_a[204], pass_b[216], pass_c[12][6], pass_s[12][6], final_R[229]);
	PE2(4, pass_a[183], pass_b[194], pass_c[11][7], pass_s[11][7], pass_a[195], pass_b[206], pass_c[11][8], pass_s[11][8], final_R[218]);
	PE2(3, pass_a[172], pass_b[182], pass_c[10][9], pass_s[10][9], pass_a[185], pass_b[195], pass_c[10][10], pass_s[10][10], final_R[206]);
	PE2(2, pass_a[160], pass_b[169], pass_c[9][11], pass_s[9][11], pass_a[174], pass_b[183], pass_c[9][12], pass_s[9][12], final_R[193]);
	PE2_tail(pass_a[147], pass_b[155], pass_c[8][13], pass_s[8][13], pass_b[170], final_R[179]);

	PE1(8, out_a[15], pass_b[240], out_a[16], PE1_pass_c[16], PE1_pass_s[16], out_c_Q[16], out_s_Q[16], final_R[264]);
	PE2(7, pass_a[218], pass_b[233], pass_c[15][0], pass_s[15][0], pass_a[226], pass_b[241], pass_c[15][1], pass_s[15][1], final_R[257]);
	PE2(6, pass_a[211], pass_b[225], pass_c[14][2], pass_s[14][2], pass_a[220], pass_b[234], pass_c[14][3], pass_s[14][3], final_R[249]);
	PE2(5, pass_a[203], pass_b[216], pass_c[13][4], pass_s[13][4], pass_a[213], pass_b[226], pass_c[13][5], pass_s[13][5], final_R[240]);
	PE2(4, pass_a[194], pass_b[206], pass_c[12][6], pass_s[12][6], pass_a[205], pass_b[217], pass_c[12][7], pass_s[12][7], final_R[230]);
	PE2(3, pass_a[184], pass_b[195], pass_c[11][8], pass_s[11][8], pass_a[196], pass_b[207], pass_c[11][9], pass_s[11][9], final_R[219]);
	PE2(2, pass_a[173], pass_b[183], pass_c[10][10], pass_s[10][10], pass_a[186], pass_b[196], pass_c[10][11], pass_s[10][11], final_R[207]);
	PE2_tail(pass_a[161], pass_b[170], pass_c[9][12], pass_s[9][12], pass_b[184], final_R[194]);

	PE3(16,out_c_Q[16], out_s_Q[16],Q_left[15], Q_right[15],Q_left[16], Q_right[16], final_Q_left[16], final_Q_right[16]);
	PE2(7, pass_a[225], pass_b[241], PE1_pass_c[16], PE1_pass_s[16], pass_a[232], pass_b[248], pass_c[16][0], pass_s[16][0], final_R[265]);
	PE2(6, pass_a[219], pass_b[234], pass_c[15][1], pass_s[15][1], pass_a[227], pass_b[242], pass_c[15][2], pass_s[15][2], final_R[258]);
	PE2(5, pass_a[212], pass_b[226], pass_c[14][3], pass_s[14][3], pass_a[221], pass_b[235], pass_c[14][4], pass_s[14][4], final_R[250]);
	PE2(4, pass_a[204], pass_b[217], pass_c[13][5], pass_s[13][5], pass_a[214], pass_b[227], pass_c[13][6], pass_s[13][6], final_R[241]);
	PE2(3, pass_a[195], pass_b[207], pass_c[12][7], pass_s[12][7], pass_a[206], pass_b[218], pass_c[12][8], pass_s[12][8], final_R[231]);
	PE2(2, pass_a[185], pass_b[196], pass_c[11][9], pass_s[11][9], pass_a[197], pass_b[208], pass_c[11][10], pass_s[11][10], final_R[220]);
	PE2_tail(pass_a[174], pass_b[184], pass_c[10][11], pass_s[10][11], pass_b[197], final_R[208]);

	PE1(7, out_a[16], pass_b[248], out_a[17], PE1_pass_c[17], PE1_pass_s[17], out_c_Q[17], out_s_Q[17], final_R[272]);
	PE2(6, pass_a[226], pass_b[242], pass_c[16][0], pass_s[16][0], pass_a[233], pass_b[249], pass_c[16][1], pass_s[16][1], final_R[266]);
	PE2(5, pass_a[220], pass_b[235], pass_c[15][2], pass_s[15][2], pass_a[228], pass_b[243], pass_c[15][3], pass_s[15][3], final_R[259]);
	PE2(4, pass_a[213], pass_b[227], pass_c[14][4], pass_s[14][4], pass_a[222], pass_b[236], pass_c[14][5], pass_s[14][5], final_R[251]);
	PE2(3, pass_a[205], pass_b[218], pass_c[13][6], pass_s[13][6], pass_a[215], pass_b[228], pass_c[13][7], pass_s[13][7], final_R[242]);
	PE2(2, pass_a[196], pass_b[208], pass_c[12][8], pass_s[12][8], pass_a[207], pass_b[219], pass_c[12][9], pass_s[12][9], final_R[232]);
	PE2_tail(pass_a[186], pass_b[197], pass_c[11][10], pass_s[11][10], pass_b[209], final_R[221]);

	PE3(17,out_c_Q[17], out_s_Q[17],Q_left[16], Q_right[16],Q_left[17], Q_right[17], final_Q_left[17], final_Q_right[17]);
	PE2(6, pass_a[232], pass_b[249], PE1_pass_c[17], PE1_pass_s[17], pass_a[238], pass_b[255], pass_c[17][0], pass_s[17][0], final_R[273]);
	PE2(5, pass_a[227], pass_b[243], pass_c[16][1], pass_s[16][1], pass_a[234], pass_b[250], pass_c[16][2], pass_s[16][2], final_R[267]);
	PE2(4, pass_a[221], pass_b[236], pass_c[15][3], pass_s[15][3], pass_a[229], pass_b[244], pass_c[15][4], pass_s[15][4], final_R[260]);
	PE2(3, pass_a[214], pass_b[228], pass_c[14][5], pass_s[14][5], pass_a[223], pass_b[237], pass_c[14][6], pass_s[14][6], final_R[252]);
	PE2(2, pass_a[206], pass_b[219], pass_c[13][7], pass_s[13][7], pass_a[216], pass_b[229], pass_c[13][8], pass_s[13][8], final_R[243]);
	PE2_tail(pass_a[197], pass_b[209], pass_c[12][9], pass_s[12][9], pass_b[220], final_R[233]);

	PE1(6, out_a[17], pass_b[255], out_a[18], PE1_pass_c[18], PE1_pass_s[18], out_c_Q[18], out_s_Q[18], final_R[279]);
	PE2(5, pass_a[233], pass_b[250], pass_c[17][0], pass_s[17][0], pass_a[239], pass_b[256], pass_c[17][1], pass_s[17][1], final_R[274]);
	PE2(4, pass_a[228], pass_b[244], pass_c[16][2], pass_s[16][2], pass_a[235], pass_b[251], pass_c[16][3], pass_s[16][3], final_R[268]);
	PE2(3, pass_a[222], pass_b[237], pass_c[15][4], pass_s[15][4], pass_a[230], pass_b[245], pass_c[15][5], pass_s[15][5], final_R[261]);
	PE2(2, pass_a[215], pass_b[229], pass_c[14][6], pass_s[14][6], pass_a[224], pass_b[238], pass_c[14][7], pass_s[14][7], final_R[253]);
	PE2_tail(pass_a[207], pass_b[220], pass_c[13][8], pass_s[13][8], pass_b[230], final_R[244]);

	PE3(18,out_c_Q[18], out_s_Q[18],Q_left[17], Q_right[17],Q_left[18], Q_right[18], final_Q_left[18], final_Q_right[18]);
	PE2(5, pass_a[238], pass_b[256], PE1_pass_c[18], PE1_pass_s[18], pass_a[243], pass_b[261], pass_c[18][0], pass_s[18][0], final_R[280]);
	PE2(4, pass_a[234], pass_b[251], pass_c[17][1], pass_s[17][1], pass_a[240], pass_b[257], pass_c[17][2], pass_s[17][2], final_R[275]);
	PE2(3, pass_a[229], pass_b[245], pass_c[16][3], pass_s[16][3], pass_a[236], pass_b[252], pass_c[16][4], pass_s[16][4], final_R[269]);
	PE2(2, pass_a[223], pass_b[238], pass_c[15][5], pass_s[15][5], pass_a[231], pass_b[246], pass_c[15][6], pass_s[15][6], final_R[262]);
	PE2_tail(pass_a[216], pass_b[230], pass_c[14][7], pass_s[14][7], pass_b[239], final_R[254]);

	PE1(5, out_a[18], pass_b[261], out_a[19], PE1_pass_c[19], PE1_pass_s[19], out_c_Q[19], out_s_Q[19], final_R[285]);
	PE2(4, pass_a[239], pass_b[257], pass_c[18][0], pass_s[18][0], pass_a[244], pass_b[262], pass_c[18][1], pass_s[18][1], final_R[281]);
	PE2(3, pass_a[235], pass_b[252], pass_c[17][2], pass_s[17][2], pass_a[241], pass_b[258], pass_c[17][3], pass_s[17][3], final_R[276]);
	PE2(2, pass_a[230], pass_b[246], pass_c[16][4], pass_s[16][4], pass_a[237], pass_b[253], pass_c[16][5], pass_s[16][5], final_R[270]);
	PE2_tail(pass_a[224], pass_b[239], pass_c[15][6], pass_s[15][6], pass_b[247], final_R[263]);

	PE3(19,out_c_Q[19], out_s_Q[19],Q_left[18], Q_right[18],Q_left[19], Q_right[19], final_Q_left[19], final_Q_right[19]);
	PE2(4, pass_a[243], pass_b[262], PE1_pass_c[19], PE1_pass_s[19], pass_a[247], pass_b[266], pass_c[19][0], pass_s[19][0], final_R[286]);
	PE2(3, pass_a[240], pass_b[258], pass_c[18][1], pass_s[18][1], pass_a[245], pass_b[263], pass_c[18][2], pass_s[18][2], final_R[282]);
	PE2(2, pass_a[236], pass_b[253], pass_c[17][3], pass_s[17][3], pass_a[242], pass_b[259], pass_c[17][4], pass_s[17][4], final_R[277]);
	PE2_tail(pass_a[231], pass_b[247], pass_c[16][5], pass_s[16][5], pass_b[254], final_R[271]);

	PE1(4, out_a[19], pass_b[266], out_a[20], PE1_pass_c[20], PE1_pass_s[20], out_c_Q[20], out_s_Q[20], final_R[290]);
	PE2(3, pass_a[244], pass_b[263], pass_c[19][0], pass_s[19][0], pass_a[248], pass_b[267], pass_c[19][1], pass_s[19][1], final_R[287]);
	PE2(2, pass_a[241], pass_b[259], pass_c[18][2], pass_s[18][2], pass_a[246], pass_b[264], pass_c[18][3], pass_s[18][3], final_R[283]);
	PE2_tail(pass_a[237], pass_b[254], pass_c[17][4], pass_s[17][4], pass_b[260], final_R[278]);

	PE3(20,out_c_Q[20], out_s_Q[20],Q_left[19], Q_right[19],Q_left[20], Q_right[20], final_Q_left[20], final_Q_right[20]);
	PE2(3, pass_a[247], pass_b[267], PE1_pass_c[20], PE1_pass_s[20], pass_a[250], pass_b[270], pass_c[20][0], pass_s[20][0], final_R[291]);
	PE2(2, pass_a[245], pass_b[264], pass_c[19][1], pass_s[19][1], pass_a[249], pass_b[268], pass_c[19][2], pass_s[19][2], final_R[288]);
	PE2_tail(pass_a[242], pass_b[260], pass_c[18][3], pass_s[18][3], pass_b[265], final_R[284]);

	PE1(3, out_a[20], pass_b[270], out_a[21], PE1_pass_c[21], PE1_pass_s[21], out_c_Q[21], out_s_Q[21], final_R[294]);
	PE2(2, pass_a[248], pass_b[268], pass_c[20][0], pass_s[20][0], pass_a[251], pass_b[271], pass_c[20][1], pass_s[20][1], final_R[292]);
	PE2_tail(pass_a[246], pass_b[265], pass_c[19][2], pass_s[19][2], pass_b[269], final_R[289]);

	PE3(21,out_c_Q[21], out_s_Q[21],Q_left[20], Q_right[20],Q_left[21], Q_right[21], final_Q_left[21], final_Q_right[21]);
	PE2(2, pass_a[250], pass_b[271], PE1_pass_c[21], PE1_pass_s[21], pass_a[252], pass_b[273], pass_c[21][0], pass_s[21][0], final_R[295]);
	PE2_tail(pass_a[249], pass_b[269], pass_c[20][1], pass_s[20][1], pass_b[272], final_R[293]);

	PE1(2, out_a[21], pass_b[273], out_a[22], PE1_pass_c[22], PE1_pass_s[22], out_c_Q[22], out_s_Q[22], final_R[297]);
	PE2_tail(pass_a[251], pass_b[272], pass_c[21][0], pass_s[21][0], pass_b[274], final_R[296]);

	PE3(22,out_c_Q[22], out_s_Q[22],Q_left[21], Q_right[21],Q_left[22], Q_right[22], final_Q_left[22], final_Q_right[22]);
	PE2_tail(pass_a[252], pass_b[274], PE1_pass_c[22], PE1_pass_s[22], pass_b[275], final_R[298]);

	PE1_tail(out_a[22], pass_b[275], out_c_Q[23], out_s_Q[23], final_R[299]);

	PE3_tail(23,out_c_Q[23], out_s_Q[23], Q_left[22], Q_right[22], final_Q_left[23], final_Q_right[23]);

// REPLACE END 
}


int top(
	stream<MATRIX_T>&A1,
	stream<MATRIX_T>&A2,
	stream<fixed_cs>&Q_L,
	stream<fixed_cs>&Q_R,
	stream<MATRIX_T>&R
){
#pragma HLS dataflow
#pragma HLS INTERFACE mode=ap_ctrl_hs port=return

#pragma HLS INTERFACE mode=axis port=A1
#pragma HLS INTERFACE mode=axis port=A2
#pragma HLS INTERFACE mode=axis port=Q_L
#pragma HLS INTERFACE mode=axis port=Q_R
#pragma HLS INTERFACE mode=axis port=R
#pragma HLS INTERFACE mode=s_axilite port=return
	stream<MATRIX_T,COLS> Feed_A[COLS] ,Feed_B[COLS];
	stream<fixed_cs,COLS*2> final_Q_left[COLS], final_Q_right[COLS];
	stream<MATRIX_T,2> final_R[LEN];

	feeder(A1, A2, Feed_A, Feed_B);
	region(Feed_A,Feed_B,final_Q_left,final_Q_right,final_R);
	collector_Q(final_Q_left, final_Q_right, Q_L, Q_R);
	collector_R(final_R, R);
	return 0;
}
  
