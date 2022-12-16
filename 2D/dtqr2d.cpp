#include "dtqr2d.h"

#include <iostream>

const int ROWS_C = ROWS;
const int COLS_C = COLS;
const int COLS_C_N = COLS_C - 1;
const int LEN_S_C = LEN_S;
const fixed_cs_mag one = 1;
using namespace std;

struct cs_bag{
	fixed_cs c , s;
};
 
MATRIX_T qrf_mag_B(MATRIX_T a, MATRIX_T b, cs_bag &cs)
{
#pragma HLS INLINE
#pragma HLS aggregate variable=cs
	fixed_double temp = a*a + b*b;;
	fixed_u mag = hls::sqrt(temp);
	fixed_cs_mag rmag = one/mag;
	cs.c = a*rmag; 
	cs.s = b*rmag;
	return mag;
}

void qrf_mm_Q(cs_bag cs, fixed_cs &op1, fixed_cs &op2)
{
#pragma HLS inline
#pragma HLS aggregate variable=cs
	fixed_cs a = op2 * cs.s + op1 * cs.c; 
	fixed_cs b = op2 * cs.c - op1 * cs.s;
	op1 = a;
	op2 = b;
}

void qrf_mm(cs_bag cs, MATRIX_T &op1, MATRIX_T &op2)
{
#pragma HLS inline
#pragma HLS aggregate variable=cs
	MATRIX_T a = op2 * cs.s + op1 * cs.c;
	MATRIX_T b = op2 * cs.c - op1 * cs.s;
	op1 = a;
	op2 = b;
}


void feeder(
	hls::stream<MATRIX_T> &A1,
	hls::stream<MATRIX_T> &A2,
	hls::stream<MATRIX_T> *Feed_A,
	hls::stream<MATRIX_T> *Feed_B)
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
	hls::stream<MATRIX_T> &in_a,
	hls::stream<MATRIX_T> &in_b,
	hls::stream<MATRIX_T> &out_a,
	hls::stream<cs_bag> & out_cs,
	hls::stream<cs_bag> & out_cs_Q,
	hls::stream<MATRIX_T> &out_R)
{
#pragma HLS FUNCTION_INSTANTIATE variable=output_size
cout<<"PE1 "<<COLS-output_size<<" start" <<endl;
	MATRIX_T A, B;
	cs_bag cs[COLS]; 
	MATRIX_T mag;
	int p = 0;
loading_Loop:
	for (int i = 0; i < COLS_C; i++)
	{
#pragma HLS loop_tripcount min=1 max=COLS_C avg=COLS_C/2
#pragma HLS pipeline II = 1
		A = in_a.read();
		B = in_b.read();
		mag = qrf_mag_B(A, B,cs[i]);
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
		out_cs_Q.write(cs[i]);
		if (i != output_size-1){
			out_cs.write(cs[i]);
		}else{
			break;
		}
	}
}

void PE1_tail(
	hls::stream<MATRIX_T> &in_a,
	hls::stream<MATRIX_T> &in_b,
	hls::stream<cs_bag> &out_cs_Q,
	hls::stream<MATRIX_T> &out_R)
{
cout<<"PE1_tail start" <<endl;
	MATRIX_T A = in_a.read();
	MATRIX_T B = in_b.read();
	cs_bag cs; 
	MATRIX_T mag = qrf_mag_B(A, B, cs);
	out_cs_Q.write(cs);
	out_R.write(mag);
}

void PE2(
	medium_int output_size, //output size = COLS - idx - idy - 1
	hls::stream<MATRIX_T> &in_a,
	hls::stream<MATRIX_T> &in_b,
	hls::stream<cs_bag> &in_cs,

	hls::stream<MATRIX_T> &pass_a,
	hls::stream<MATRIX_T> &pass_b,
	hls::stream<cs_bag> &pass_cs,

	hls::stream<MATRIX_T> &final_R)

{
#pragma HLS FUNCTION_INSTANTIATE variable=output_size
cout<<"PE2 "<<COLS-output_size<<" start" <<endl;
	MATRIX_T A, B;
	cs_bag cs;
	medium_int pass_cs_range = output_size-1;
	for (medium_int i = 0; i < COLS_C; i++)
	{
// #pragma HLS pipeline II = 1
		if (i < output_size){
			int a = in_cs.empty();
			cs = in_cs.read();
			if (i != pass_cs_range){
				pass_cs.write(cs);
			}
			A = in_a.read();
			B= in_b.read();
			qrf_mm(cs, A, B);
			if (i == 0)
				final_R.write(A);
			else
				pass_a.write(A);
			pass_b.write(B);
		}
	}
	cout<<"PE2 "<<COLS-output_size<<" end" <<endl;

}

void PE2_tail(
	hls::stream<MATRIX_T> &in_a,
	hls::stream<MATRIX_T> &in_b,
	hls::stream<cs_bag> &in_cs,
	hls::stream<MATRIX_T> &out_b,
	hls::stream<MATRIX_T> &out_a_R)
{
cout<<"PE2 tail start" <<endl;
	MATRIX_T A, B;
	cs_bag cs;
	cs = in_cs.read();
	A = in_a.read();
	B = in_b.read();
	qrf_mm(cs, A, B);
	out_a_R.write(A);
	out_b.write(B);
}

void PE3_head(hls::stream<cs_bag> &in_cs,
			  hls::stream<fixed_cs> &out_Q_left,
			  hls::stream<fixed_cs> &out_Q_right,
			  hls::stream<fixed_cs> &final_Q_left,
			  hls::stream<fixed_cs> &final_Q_right)
{
cout<<"PE3 head start" <<endl;
	cs_bag cs;
	fixed_cs x, y;
	for (medium_int i = 0; i < COLS_C; i++)
	{
		cs = in_cs.read();
		for (small_int j = 0; j < 2; j++)
		{
			if (j==0){
				x = 1;
				y = 0;
			}else{
				x = 0;
				y = 1;
			}
			qrf_mm_Q(cs, x, y);
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
		 hls::stream<cs_bag> &in_cs,
		 hls::stream<fixed_cs> &in_Q_left,
		 hls::stream<fixed_cs> &in_Q_right,
		 hls::stream<fixed_cs> &out_Q_left,
		 hls::stream<fixed_cs> &out_Q_right,
		 hls::stream<fixed_cs> &final_Q_left,
		 hls::stream<fixed_cs> &final_Q_right)
{
cout<<"PE3 "<<id<<" start" <<endl;

#pragma HLS FUNCTION_INSTANTIATE variable=id
	cs_bag cs;
	fixed_cs x, y;
	for (medium_int i = 0; i < COLS; i++)
	{
		if (i < COLS-id){
		cs = in_cs.read();
		for (small_int _ = 0; _ < 2; _++)
		{
			for (medium_int j = 0; j < COLS; j++)
			{
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

					qrf_mm_Q(cs, x, y);

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
	hls::stream<cs_bag> &in_cs,
	hls::stream<fixed_cs> &in_Q_left,
	hls::stream<fixed_cs> &in_Q_right,
	hls::stream<fixed_cs> &final_Q_left,
	hls::stream<fixed_cs> &final_Q_right)
{
	cout<<"PE3_tail "<<id<<" start" <<endl;
	cs_bag cs;
	fixed_cs x, y;
	cs = in_cs.read();
	for (small_int _ = 0; _ < 2; _++)
	{
		x = 0;
		y = in_Q_right.read();
		qrf_mm_Q(cs, x, y);
		final_Q_left.write(x);
		final_Q_right.write(y);
		for (medium_int j = 0; j < COLS; j++)
		{
			if (j < id-1)
			{
				x = in_Q_left.read();
				y = in_Q_right.read();
				qrf_mm_Q(cs, x, y);
				final_Q_left.write(x);
				final_Q_right.write(y);
			}
		}
		x = in_Q_left.read();
		y = 0;
		qrf_mm_Q(cs, x, y);
		final_Q_left.write(x);
		final_Q_right.write(y);
	}
}

void collector_Q(
	hls::stream<fixed_cs> *final_Q_left,
	hls::stream<fixed_cs> *final_Q_right,
	hls::stream<fixed_cs> &Q_L,
	hls::stream<fixed_cs> &Q_R)
{
	cout<<"collQ start" <<endl;
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
	cout<<"collQ end" <<endl;

}

void collector_R(hls::stream<MATRIX_T> *in_R, hls::stream<MATRIX_T> &R)
{
	cout<<"collR start" <<endl;
	for (int i = 0; i < LEN; i++)
	{
		R.write(in_R[i].read());
	}
	cout<<"collR end" <<endl;

}

void region(
	hls::stream<MATRIX_T,COLS> *Feed_A,
	hls::stream<MATRIX_T,COLS> *Feed_B,
	hls::stream<fixed_cs,COLS*2> *final_Q_left,
	hls::stream<fixed_cs,COLS*2> *final_Q_right,
	hls::stream<MATRIX_T,2> *final_R
){
	#pragma HLS dataflow
	#pragma HLS INTERFACE mode=ap_ctrl_none port=return

	hls::stream<MATRIX_T,COLS> out_a[COLS]; 
	hls::stream<cs_bag,COLS> PE1_pass_cs[COLS];
	hls::stream<cs_bag,2> pass_cs[COLS][COLS];  
	hls::stream<MATRIX_T,COLS> pass_a[LEN_S], pass_b[LEN_S];
	hls::stream<cs_bag,COLS> out_cs_Q[COLS];
	hls::stream<fixed_cs,COLS*2> Q_left[COLS], Q_right[COLS]; 

// REPLACE START
	PE1(12, Feed_A[0], Feed_B[0], out_a[0], PE1_pass_cs[0], out_cs_Q[0], final_R[0]);

	PE3_head(out_cs_Q[0], Q_left[0], Q_right[0], final_Q_left[0], final_Q_right[0]);
	PE2(11, Feed_A[1], Feed_B[1], PE1_pass_cs[0], pass_a[0], pass_b[0], pass_cs[0][0], final_R[1]);

	PE1(11, out_a[0], pass_b[0], out_a[1], PE1_pass_cs[1], out_cs_Q[1], final_R[12]);
	PE2(10, Feed_A[2], Feed_B[2], pass_cs[0][0], pass_a[1], pass_b[1], pass_cs[0][1], final_R[2]);

	PE3(1,out_cs_Q[1],Q_left[0], Q_right[0],Q_left[1], Q_right[1], final_Q_left[1], final_Q_right[1]);
	PE2(10, pass_a[0], pass_b[1], PE1_pass_cs[1], pass_a[10], pass_b[11], pass_cs[1][0], final_R[13]);
	PE2(9, Feed_A[3], Feed_B[3], pass_cs[0][1], pass_a[2], pass_b[2], pass_cs[0][2], final_R[3]);

	PE1(10, out_a[1], pass_b[11], out_a[2], PE1_pass_cs[2], out_cs_Q[2], final_R[23]);
	PE2(9, pass_a[1], pass_b[2], pass_cs[1][0], pass_a[11], pass_b[12], pass_cs[1][1], final_R[14]);
	PE2(8, Feed_A[4], Feed_B[4], pass_cs[0][2], pass_a[3], pass_b[3], pass_cs[0][3], final_R[4]);

	PE3(2,out_cs_Q[2],Q_left[1], Q_right[1],Q_left[2], Q_right[2], final_Q_left[2], final_Q_right[2]);
	PE2(9, pass_a[10], pass_b[12], PE1_pass_cs[2], pass_a[19], pass_b[21], pass_cs[2][0], final_R[24]);
	PE2(8, pass_a[2], pass_b[3], pass_cs[1][1], pass_a[12], pass_b[13], pass_cs[1][2], final_R[15]);
	PE2(7, Feed_A[5], Feed_B[5], pass_cs[0][3], pass_a[4], pass_b[4], pass_cs[0][4], final_R[5]);

	PE1(9, out_a[2], pass_b[21], out_a[3], PE1_pass_cs[3], out_cs_Q[3], final_R[33]);
	PE2(8, pass_a[11], pass_b[13], pass_cs[2][0], pass_a[20], pass_b[22], pass_cs[2][1], final_R[25]);
	PE2(7, pass_a[3], pass_b[4], pass_cs[1][2], pass_a[13], pass_b[14], pass_cs[1][3], final_R[16]);
	PE2(6, Feed_A[6], Feed_B[6], pass_cs[0][4], pass_a[5], pass_b[5], pass_cs[0][5], final_R[6]);

	PE3(3,out_cs_Q[3],Q_left[2], Q_right[2],Q_left[3], Q_right[3], final_Q_left[3], final_Q_right[3]);
	PE2(8, pass_a[19], pass_b[22], PE1_pass_cs[3], pass_a[27], pass_b[30], pass_cs[3][0], final_R[34]);
	PE2(7, pass_a[12], pass_b[14], pass_cs[2][1], pass_a[21], pass_b[23], pass_cs[2][2], final_R[26]);
	PE2(6, pass_a[4], pass_b[5], pass_cs[1][3], pass_a[14], pass_b[15], pass_cs[1][4], final_R[17]);
	PE2(5, Feed_A[7], Feed_B[7], pass_cs[0][5], pass_a[6], pass_b[6], pass_cs[0][6], final_R[7]);

	PE1(8, out_a[3], pass_b[30], out_a[4], PE1_pass_cs[4], out_cs_Q[4], final_R[42]);
	PE2(7, pass_a[20], pass_b[23], pass_cs[3][0], pass_a[28], pass_b[31], pass_cs[3][1], final_R[35]);
	PE2(6, pass_a[13], pass_b[15], pass_cs[2][2], pass_a[22], pass_b[24], pass_cs[2][3], final_R[27]);
	PE2(5, pass_a[5], pass_b[6], pass_cs[1][4], pass_a[15], pass_b[16], pass_cs[1][5], final_R[18]);
	PE2(4, Feed_A[8], Feed_B[8], pass_cs[0][6], pass_a[7], pass_b[7], pass_cs[0][7], final_R[8]);

	PE3(4,out_cs_Q[4],Q_left[3], Q_right[3],Q_left[4], Q_right[4], final_Q_left[4], final_Q_right[4]);
	PE2(7, pass_a[27], pass_b[31], PE1_pass_cs[4], pass_a[34], pass_b[38], pass_cs[4][0], final_R[43]);
	PE2(6, pass_a[21], pass_b[24], pass_cs[3][1], pass_a[29], pass_b[32], pass_cs[3][2], final_R[36]);
	PE2(5, pass_a[14], pass_b[16], pass_cs[2][3], pass_a[23], pass_b[25], pass_cs[2][4], final_R[28]);
	PE2(4, pass_a[6], pass_b[7], pass_cs[1][5], pass_a[16], pass_b[17], pass_cs[1][6], final_R[19]);
	PE2(3, Feed_A[9], Feed_B[9], pass_cs[0][7], pass_a[8], pass_b[8], pass_cs[0][8], final_R[9]);

	PE1(7, out_a[4], pass_b[38], out_a[5], PE1_pass_cs[5], out_cs_Q[5], final_R[50]);
	PE2(6, pass_a[28], pass_b[32], pass_cs[4][0], pass_a[35], pass_b[39], pass_cs[4][1], final_R[44]);
	PE2(5, pass_a[22], pass_b[25], pass_cs[3][2], pass_a[30], pass_b[33], pass_cs[3][3], final_R[37]);
	PE2(4, pass_a[15], pass_b[17], pass_cs[2][4], pass_a[24], pass_b[26], pass_cs[2][5], final_R[29]);
	PE2(3, pass_a[7], pass_b[8], pass_cs[1][6], pass_a[17], pass_b[18], pass_cs[1][7], final_R[20]);
	PE2(2, Feed_A[10], Feed_B[10], pass_cs[0][8], pass_a[9], pass_b[9], pass_cs[0][9], final_R[10]);

	PE3(5,out_cs_Q[5],Q_left[4], Q_right[4],Q_left[5], Q_right[5], final_Q_left[5], final_Q_right[5]);
	PE2(6, pass_a[34], pass_b[39], PE1_pass_cs[5], pass_a[40], pass_b[45], pass_cs[5][0], final_R[51]);
	PE2(5, pass_a[29], pass_b[33], pass_cs[4][1], pass_a[36], pass_b[40], pass_cs[4][2], final_R[45]);
	PE2(4, pass_a[23], pass_b[26], pass_cs[3][3], pass_a[31], pass_b[34], pass_cs[3][4], final_R[38]);
	PE2(3, pass_a[16], pass_b[18], pass_cs[2][5], pass_a[25], pass_b[27], pass_cs[2][6], final_R[30]);
	PE2(2, pass_a[8], pass_b[9], pass_cs[1][7], pass_a[18], pass_b[19], pass_cs[1][8], final_R[21]);
	PE2_tail(Feed_A[11], Feed_B[11], pass_cs[0][9], pass_b[10], final_R[11]);

	PE1(6, out_a[5], pass_b[45], out_a[6], PE1_pass_cs[6], out_cs_Q[6], final_R[57]);
	PE2(5, pass_a[35], pass_b[40], pass_cs[5][0], pass_a[41], pass_b[46], pass_cs[5][1], final_R[52]);
	PE2(4, pass_a[30], pass_b[34], pass_cs[4][2], pass_a[37], pass_b[41], pass_cs[4][3], final_R[46]);
	PE2(3, pass_a[24], pass_b[27], pass_cs[3][4], pass_a[32], pass_b[35], pass_cs[3][5], final_R[39]);
	PE2(2, pass_a[17], pass_b[19], pass_cs[2][6], pass_a[26], pass_b[28], pass_cs[2][7], final_R[31]);
	PE2_tail(pass_a[9], pass_b[10], pass_cs[1][8], pass_b[20], final_R[22]);

	PE3(6,out_cs_Q[6],Q_left[5], Q_right[5],Q_left[6], Q_right[6], final_Q_left[6], final_Q_right[6]);
	PE2(5, pass_a[40], pass_b[46], PE1_pass_cs[6], pass_a[45], pass_b[51], pass_cs[6][0], final_R[58]);
	PE2(4, pass_a[36], pass_b[41], pass_cs[5][1], pass_a[42], pass_b[47], pass_cs[5][2], final_R[53]);
	PE2(3, pass_a[31], pass_b[35], pass_cs[4][3], pass_a[38], pass_b[42], pass_cs[4][4], final_R[47]);
	PE2(2, pass_a[25], pass_b[28], pass_cs[3][5], pass_a[33], pass_b[36], pass_cs[3][6], final_R[40]);
	PE2_tail(pass_a[18], pass_b[20], pass_cs[2][7], pass_b[29], final_R[32]);

	PE1(5, out_a[6], pass_b[51], out_a[7], PE1_pass_cs[7], out_cs_Q[7], final_R[63]);
	PE2(4, pass_a[41], pass_b[47], pass_cs[6][0], pass_a[46], pass_b[52], pass_cs[6][1], final_R[59]);
	PE2(3, pass_a[37], pass_b[42], pass_cs[5][2], pass_a[43], pass_b[48], pass_cs[5][3], final_R[54]);
	PE2(2, pass_a[32], pass_b[36], pass_cs[4][4], pass_a[39], pass_b[43], pass_cs[4][5], final_R[48]);
	PE2_tail(pass_a[26], pass_b[29], pass_cs[3][6], pass_b[37], final_R[41]);

	PE3(7,out_cs_Q[7],Q_left[6], Q_right[6],Q_left[7], Q_right[7], final_Q_left[7], final_Q_right[7]);
	PE2(4, pass_a[45], pass_b[52], PE1_pass_cs[7], pass_a[49], pass_b[56], pass_cs[7][0], final_R[64]);
	PE2(3, pass_a[42], pass_b[48], pass_cs[6][1], pass_a[47], pass_b[53], pass_cs[6][2], final_R[60]);
	PE2(2, pass_a[38], pass_b[43], pass_cs[5][3], pass_a[44], pass_b[49], pass_cs[5][4], final_R[55]);
	PE2_tail(pass_a[33], pass_b[37], pass_cs[4][5], pass_b[44], final_R[49]);

	PE1(4, out_a[7], pass_b[56], out_a[8], PE1_pass_cs[8], out_cs_Q[8], final_R[68]);
	PE2(3, pass_a[46], pass_b[53], pass_cs[7][0], pass_a[50], pass_b[57], pass_cs[7][1], final_R[65]);
	PE2(2, pass_a[43], pass_b[49], pass_cs[6][2], pass_a[48], pass_b[54], pass_cs[6][3], final_R[61]);
	PE2_tail(pass_a[39], pass_b[44], pass_cs[5][4], pass_b[50], final_R[56]);

	PE3(8,out_cs_Q[8],Q_left[7], Q_right[7],Q_left[8], Q_right[8], final_Q_left[8], final_Q_right[8]);
	PE2(3, pass_a[49], pass_b[57], PE1_pass_cs[8], pass_a[52], pass_b[60], pass_cs[8][0], final_R[69]);
	PE2(2, pass_a[47], pass_b[54], pass_cs[7][1], pass_a[51], pass_b[58], pass_cs[7][2], final_R[66]);
	PE2_tail(pass_a[44], pass_b[50], pass_cs[6][3], pass_b[55], final_R[62]);

	PE1(3, out_a[8], pass_b[60], out_a[9], PE1_pass_cs[9], out_cs_Q[9], final_R[72]);
	PE2(2, pass_a[50], pass_b[58], pass_cs[8][0], pass_a[53], pass_b[61], pass_cs[8][1], final_R[70]);
	PE2_tail(pass_a[48], pass_b[55], pass_cs[7][2], pass_b[59], final_R[67]);

	PE3(9,out_cs_Q[9],Q_left[8], Q_right[8],Q_left[9], Q_right[9], final_Q_left[9], final_Q_right[9]);
	PE2(2, pass_a[52], pass_b[61], PE1_pass_cs[9], pass_a[54], pass_b[63], pass_cs[9][0], final_R[73]);
	PE2_tail(pass_a[51], pass_b[59], pass_cs[8][1], pass_b[62], final_R[71]);

	PE1(2, out_a[9], pass_b[63], out_a[10], PE1_pass_cs[10], out_cs_Q[10], final_R[75]);
	PE2_tail(pass_a[53], pass_b[62], pass_cs[9][0], pass_b[64], final_R[74]);

	PE3(10,out_cs_Q[10],Q_left[9], Q_right[9],Q_left[10], Q_right[10], final_Q_left[10], final_Q_right[10]);
	PE2_tail(pass_a[54], pass_b[64], PE1_pass_cs[10], pass_b[65], final_R[76]);

	PE1_tail(out_a[10], pass_b[65], out_cs_Q[11], final_R[77]);

	PE3_tail(11,out_cs_Q[11], Q_left[10], Q_right[10], final_Q_left[11], final_Q_right[11]);

// REPLACE END 
}


int top(
	hls::stream<MATRIX_T>&A1,
	hls::stream<MATRIX_T>&A2,
	hls::stream<fixed_cs>&Q_L,
	hls::stream<fixed_cs>&Q_R,
	hls::stream<MATRIX_T>&R
){
#pragma HLS dataflow
#pragma HLS INTERFACE mode=ap_ctrl_hs port=return

#pragma HLS INTERFACE mode=axis port=A1
#pragma HLS INTERFACE mode=axis port=A2
#pragma HLS INTERFACE mode=axis port=Q_L
#pragma HLS INTERFACE mode=axis port=Q_R
#pragma HLS INTERFACE mode=axis port=R
#pragma HLS INTERFACE mode=s_axilite port=return

	hls::stream<MATRIX_T,COLS> Feed_A[COLS] ,Feed_B[COLS];
	hls::stream<fixed_cs,COLS*2> final_Q_left[COLS], final_Q_right[COLS];
	hls::stream<MATRIX_T,2> final_R[LEN];

	feeder(A1, A2, Feed_A, Feed_B);
	region(Feed_A,Feed_B,final_Q_left,final_Q_right,final_R);
	collector_Q(final_Q_left, final_Q_right, Q_L, Q_R);
	collector_R(final_R, R);
	return 0;
}
  
