#include "dtqr2d.h"

#include <iostream>

const int ROWS_C = ROWS;
const int COLS_C = COLS;
const int COLS_C_N = COLS_C - 1;
const int LEN_S_C = LEN_S;
const fixed_cs_mag one = 1;
using namespace std;
struct magreturn{
	fixed_cs c ,s;
	MATRIX_T mag;
};

struct cs_bag{
	fixed_cs c , s;
};
//
//void qrf_mag_stream(
//	hls::stream<MATRIX_T> &in_a,
//	hls::stream<MATRIX_T> &in_b,
//	hls::stream<fixed_cs> &out_c,
//	hls::stream<fixed_cs> &out_s,
//	hls::stream<MATRIX_T> &out_mag)
//{
//#pragma HLS INLINE OFF
//#pragma HLS ALLOCATION function instances=qrf_mag_hls::stream limit=3
//	MATRIX_T A = in_a.read();
//	MATRIX_T B = in_b.read();
//	fixed_double temp = A*A + B*B;;
//	fixed_u mag = hls::sqrt(temp);
//	fixed_cs_mag rmag = one/mag;
//	fixed_cs c = A*rmag;
//	fixed_cs s = B*rmag;
//	out_c.write(c);
//	out_s.write(s);
//	out_mag.write(mag);
//}
//
//magreturn qrf_mag_C(MATRIX_T a, MATRIX_T b)
//{
//#pragma HLS INLINE OFF
//	magreturn ret;
//	fixed_double temp = a*a + b*b;;
//	fixed_u mag = hls::sqrt(temp);
//	fixed_cs_mag rmag = one/mag;
//	ret.c = a*rmag;
//	ret.s = b*rmag;
//	ret.mag = mag;
//	return ret;
//}


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

	fixed_cs a = op2 * cs.s + op1 * cs.c; // multi
	fixed_cs b = op2 * cs.c - op1 * cs.s;
	op1 = a;
	op2 = b;
}

void qrf_mm(cs_bag cs, MATRIX_T &op1, MATRIX_T &op2)
{
#pragma HLS inline off
#pragma HLS aggregate variable=cs
	MATRIX_T a = op2 * cs.s + op1 * cs.c;
	MATRIX_T b = op2 * cs.c - op1 * cs.s;
	op1 = a;
	op2 = b;
}

// Feeder to prepare the data;

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
//void PE1_1(
//	medium_int output_size, // id start at 0
//	hls::stream<MATRIX_T> &in_a,
//	hls::stream<MATRIX_T> &in_b,
//	hls::stream<MATRIX_T> &out_a,
//	hls::stream<cs_bag> &out_cs,
//	hls::stream<cs_bag> &out_cs_Q,
//	hls::stream<MATRIX_T> &out_R)
//{
////#pragma HLS FUNCTION_INSTANTIATE variable=output_size
//	MATRIX_T A, B;
//	cs_bag cs[COLS];
//	MATRIX_T mag;
//	hls::stream<MATRIX_T,2> send_A[COLS] ,send_B[COLS];
//	hls::stream<fixed_cs,2> get_c[COLS], get_s[COLS];
//	hls::stream<MATRIX_T,2> get_mag[LEN];
//	int p = 0;
//loading_Loop:
//	for (int i = 0; i < COLS_C; i++)
//	{
//#pragma HLS loop_tripcount min=1 max=COLS_C avg=COLS_C/2
//#pragma HLS pipeline II = 3
//		send_A[i].write(in_a.read());
//		send_B[i].write(in_b.read());
//		qrf_mag_hls::stream(send_A[i],send_B[i],get_c[i],get_s[i],get_mag[i]);
//		c[i] = get_c[i].read();
//		s[i] = get_s[i].read();
//		if (i == 0)
//			out_R.write(get_mag[i].read());
//		else
//			out_a.write(get_mag[i].read());
//		if (i == output_size-1)
//			break;
//	}
//
//transfer_Loop_1:
//	for (int i = 0; i < COLS_C; i++)
//	{
//#pragma HLS pipeline II = 1
//		out_c_Q.write(c[i]);
//		out_s_Q.write(s[i]);
//		if (i != output_size-1){
//			out_c.write(c[i]);
//			out_s.write(s[i]);
//		}else{
//			break;
//		}
//	}
//}
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
#pragma HLS ALLOCATION function instances=qrf_mag_B limit=4
cout<<"PE1 "<<COLS-output_size<<" start" <<endl;
	MATRIX_T A, B;
	cs_bag cs[COLS]; 
	MATRIX_T mag;
	int p = 0;
loading_Loop:
	for (int i = 0; i < COLS_C; i++)
	{
#pragma HLS loop_tripcount min=1 max=COLS_C avg=COLS_C/2
#pragma HLS pipeline II = 3
		A = in_a.read();
		B = in_b.read();
		mag = qrf_mag_B(A, B,cs[i]);
//		cout<<output_size<<"  "<< cs[i].c<<endl;
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
//void PE1_2(
//	medium_int output_size, // id start at 0
//	hls::stream<MATRIX_T> &in_a,
//	hls::stream<MATRIX_T> &in_b,
//	hls::stream<MATRIX_T> &out_a,
//	hls::stream<fixed_cs> &out_c,
//	hls::stream<fixed_cs> &out_s,
//	hls::stream<fixed_cs> &out_c_Q,
//	hls::stream<fixed_cs> &out_s_Q,
//	hls::stream<MATRIX_T> &out_R)
//{
////#pragma HLS FUNCTION_INSTANTIATE variable=output_size
//#pragma HLS ALLOCATION function instances=qrf_mag_C limit=1
//
//	MATRIX_T A, B;
//	fixed_cs c[COLS], s[COLS];
//	magreturn mgr;
//	MATRIX_T mag;
//	int p = 0;
//loading_Loop:
//	for (int i = 0; i < COLS_C; i++)
//	{
//#pragma HLS loop_tripcount min=1 max=COLS_C avg=COLS_C/2
//#pragma HLS pipeline II = 3
//		A = in_a.read();
//		B = in_b.read();
//		mgr = qrf_mag_C(A, B);
//		c[i]=mgr.c;
//		s[i]=mgr.s;
//		if (i == 0)
//			out_R.write(mgr.mag);
//		else
//			out_a.write(mgr.mag);
//		if (i == output_size-1)
//			break;
//	}
//
//transfer_Loop_1:
//	for (int i = 0; i < COLS_C; i++)
//	{
//#pragma HLS pipeline II = 1
//		out_c_Q.write(c[i]);
//		out_s_Q.write(s[i]);
//		if (i != output_size-1){
//			out_c.write(c[i]);
//			out_s.write(s[i]);
//		}else{
//			break;
//		}
//	}
//}
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

cout<<"PE2 "<<COLS-output_size<<" start" <<endl;
	MATRIX_T A, B;
	cs_bag cs;
	medium_int pass_cs_range = output_size-1;
	for (medium_int i = 0; i < COLS_C; i++)
	{
#pragma HLS pipeline II = 3
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
#pragma HLS pipeline II = 1
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

// #pragma HLS FUNCTION_INSTANTIATE variable=id
	cs_bag cs;
	fixed_cs x, y;
	for (medium_int i = 0; i < COLS; i++)
	{
		if (i < COLS-id){
		cs = in_cs.read();
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
		for (medium_int j = 0; j < COLS-2; j++)
#pragma HLS pipeline II = 1
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
#pragma HLS ALLOCATION function instances=PE1 limit=4
#pragma HLS ALLOCATION function instances=PE2 limit=COLS_C


// REPLACE START
	PE1(32, Feed_A[0], Feed_B[0], out_a[0], PE1_pass_cs[0], out_cs_Q[0], final_R[0]);

	PE3_head(out_cs_Q[0], Q_left[0], Q_right[0], final_Q_left[0], final_Q_right[0]);
	PE2(31, Feed_A[1], Feed_B[1], PE1_pass_cs[0], pass_a[0], pass_b[0], pass_cs[0][0], final_R[1]);

	PE1(31, out_a[0], pass_b[0], out_a[1], PE1_pass_cs[1], out_cs_Q[1], final_R[32]);
	PE2(30, Feed_A[2], Feed_B[2], pass_cs[0][0], pass_a[1], pass_b[1], pass_cs[0][1], final_R[2]);

	PE3(1,out_cs_Q[1],Q_left[0], Q_right[0],Q_left[1], Q_right[1], final_Q_left[1], final_Q_right[1]);
	PE2(30, pass_a[0], pass_b[1], PE1_pass_cs[1], pass_a[30], pass_b[31], pass_cs[1][0], final_R[33]);
	PE2(29, Feed_A[3], Feed_B[3], pass_cs[0][1], pass_a[2], pass_b[2], pass_cs[0][2], final_R[3]);

	PE1(30, out_a[1], pass_b[31], out_a[2], PE1_pass_cs[2], out_cs_Q[2], final_R[63]);
	PE2(29, pass_a[1], pass_b[2], pass_cs[1][0], pass_a[31], pass_b[32], pass_cs[1][1], final_R[34]);
	PE2(28, Feed_A[4], Feed_B[4], pass_cs[0][2], pass_a[3], pass_b[3], pass_cs[0][3], final_R[4]);

	PE3(2,out_cs_Q[2],Q_left[1], Q_right[1],Q_left[2], Q_right[2], final_Q_left[2], final_Q_right[2]);
	PE2(29, pass_a[30], pass_b[32], PE1_pass_cs[2], pass_a[59], pass_b[61], pass_cs[2][0], final_R[64]);
	PE2(28, pass_a[2], pass_b[3], pass_cs[1][1], pass_a[32], pass_b[33], pass_cs[1][2], final_R[35]);
	PE2(27, Feed_A[5], Feed_B[5], pass_cs[0][3], pass_a[4], pass_b[4], pass_cs[0][4], final_R[5]);

	PE1(29, out_a[2], pass_b[61], out_a[3], PE1_pass_cs[3], out_cs_Q[3], final_R[93]);
	PE2(28, pass_a[31], pass_b[33], pass_cs[2][0], pass_a[60], pass_b[62], pass_cs[2][1], final_R[65]);
	PE2(27, pass_a[3], pass_b[4], pass_cs[1][2], pass_a[33], pass_b[34], pass_cs[1][3], final_R[36]);
	PE2(26, Feed_A[6], Feed_B[6], pass_cs[0][4], pass_a[5], pass_b[5], pass_cs[0][5], final_R[6]);

	PE3(3,out_cs_Q[3],Q_left[2], Q_right[2],Q_left[3], Q_right[3], final_Q_left[3], final_Q_right[3]);
	PE2(28, pass_a[59], pass_b[62], PE1_pass_cs[3], pass_a[87], pass_b[90], pass_cs[3][0], final_R[94]);
	PE2(27, pass_a[32], pass_b[34], pass_cs[2][1], pass_a[61], pass_b[63], pass_cs[2][2], final_R[66]);
	PE2(26, pass_a[4], pass_b[5], pass_cs[1][3], pass_a[34], pass_b[35], pass_cs[1][4], final_R[37]);
	PE2(25, Feed_A[7], Feed_B[7], pass_cs[0][5], pass_a[6], pass_b[6], pass_cs[0][6], final_R[7]);

	PE1(28, out_a[3], pass_b[90], out_a[4], PE1_pass_cs[4], out_cs_Q[4], final_R[122]);
	PE2(27, pass_a[60], pass_b[63], pass_cs[3][0], pass_a[88], pass_b[91], pass_cs[3][1], final_R[95]);
	PE2(26, pass_a[33], pass_b[35], pass_cs[2][2], pass_a[62], pass_b[64], pass_cs[2][3], final_R[67]);
	PE2(25, pass_a[5], pass_b[6], pass_cs[1][4], pass_a[35], pass_b[36], pass_cs[1][5], final_R[38]);
	PE2(24, Feed_A[8], Feed_B[8], pass_cs[0][6], pass_a[7], pass_b[7], pass_cs[0][7], final_R[8]);

	PE3(4,out_cs_Q[4],Q_left[3], Q_right[3],Q_left[4], Q_right[4], final_Q_left[4], final_Q_right[4]);
	PE2(27, pass_a[87], pass_b[91], PE1_pass_cs[4], pass_a[114], pass_b[118], pass_cs[4][0], final_R[123]);
	PE2(26, pass_a[61], pass_b[64], pass_cs[3][1], pass_a[89], pass_b[92], pass_cs[3][2], final_R[96]);
	PE2(25, pass_a[34], pass_b[36], pass_cs[2][3], pass_a[63], pass_b[65], pass_cs[2][4], final_R[68]);
	PE2(24, pass_a[6], pass_b[7], pass_cs[1][5], pass_a[36], pass_b[37], pass_cs[1][6], final_R[39]);
	PE2(23, Feed_A[9], Feed_B[9], pass_cs[0][7], pass_a[8], pass_b[8], pass_cs[0][8], final_R[9]);

	PE1(27, out_a[4], pass_b[118], out_a[5], PE1_pass_cs[5], out_cs_Q[5], final_R[150]);
	PE2(26, pass_a[88], pass_b[92], pass_cs[4][0], pass_a[115], pass_b[119], pass_cs[4][1], final_R[124]);
	PE2(25, pass_a[62], pass_b[65], pass_cs[3][2], pass_a[90], pass_b[93], pass_cs[3][3], final_R[97]);
	PE2(24, pass_a[35], pass_b[37], pass_cs[2][4], pass_a[64], pass_b[66], pass_cs[2][5], final_R[69]);
	PE2(23, pass_a[7], pass_b[8], pass_cs[1][6], pass_a[37], pass_b[38], pass_cs[1][7], final_R[40]);
	PE2(22, Feed_A[10], Feed_B[10], pass_cs[0][8], pass_a[9], pass_b[9], pass_cs[0][9], final_R[10]);

	PE3(5,out_cs_Q[5],Q_left[4], Q_right[4],Q_left[5], Q_right[5], final_Q_left[5], final_Q_right[5]);
	PE2(26, pass_a[114], pass_b[119], PE1_pass_cs[5], pass_a[140], pass_b[145], pass_cs[5][0], final_R[151]);
	PE2(25, pass_a[89], pass_b[93], pass_cs[4][1], pass_a[116], pass_b[120], pass_cs[4][2], final_R[125]);
	PE2(24, pass_a[63], pass_b[66], pass_cs[3][3], pass_a[91], pass_b[94], pass_cs[3][4], final_R[98]);
	PE2(23, pass_a[36], pass_b[38], pass_cs[2][5], pass_a[65], pass_b[67], pass_cs[2][6], final_R[70]);
	PE2(22, pass_a[8], pass_b[9], pass_cs[1][7], pass_a[38], pass_b[39], pass_cs[1][8], final_R[41]);
	PE2(21, Feed_A[11], Feed_B[11], pass_cs[0][9], pass_a[10], pass_b[10], pass_cs[0][10], final_R[11]);

	PE1(26, out_a[5], pass_b[145], out_a[6], PE1_pass_cs[6], out_cs_Q[6], final_R[177]);
	PE2(25, pass_a[115], pass_b[120], pass_cs[5][0], pass_a[141], pass_b[146], pass_cs[5][1], final_R[152]);
	PE2(24, pass_a[90], pass_b[94], pass_cs[4][2], pass_a[117], pass_b[121], pass_cs[4][3], final_R[126]);
	PE2(23, pass_a[64], pass_b[67], pass_cs[3][4], pass_a[92], pass_b[95], pass_cs[3][5], final_R[99]);
	PE2(22, pass_a[37], pass_b[39], pass_cs[2][6], pass_a[66], pass_b[68], pass_cs[2][7], final_R[71]);
	PE2(21, pass_a[9], pass_b[10], pass_cs[1][8], pass_a[39], pass_b[40], pass_cs[1][9], final_R[42]);
	PE2(20, Feed_A[12], Feed_B[12], pass_cs[0][10], pass_a[11], pass_b[11], pass_cs[0][11], final_R[12]);

	PE3(6,out_cs_Q[6],Q_left[5], Q_right[5],Q_left[6], Q_right[6], final_Q_left[6], final_Q_right[6]);
	PE2(25, pass_a[140], pass_b[146], PE1_pass_cs[6], pass_a[165], pass_b[171], pass_cs[6][0], final_R[178]);
	PE2(24, pass_a[116], pass_b[121], pass_cs[5][1], pass_a[142], pass_b[147], pass_cs[5][2], final_R[153]);
	PE2(23, pass_a[91], pass_b[95], pass_cs[4][3], pass_a[118], pass_b[122], pass_cs[4][4], final_R[127]);
	PE2(22, pass_a[65], pass_b[68], pass_cs[3][5], pass_a[93], pass_b[96], pass_cs[3][6], final_R[100]);
	PE2(21, pass_a[38], pass_b[40], pass_cs[2][7], pass_a[67], pass_b[69], pass_cs[2][8], final_R[72]);
	PE2(20, pass_a[10], pass_b[11], pass_cs[1][9], pass_a[40], pass_b[41], pass_cs[1][10], final_R[43]);
	PE2(19, Feed_A[13], Feed_B[13], pass_cs[0][11], pass_a[12], pass_b[12], pass_cs[0][12], final_R[13]);

	PE1(25, out_a[6], pass_b[171], out_a[7], PE1_pass_cs[7], out_cs_Q[7], final_R[203]);
	PE2(24, pass_a[141], pass_b[147], pass_cs[6][0], pass_a[166], pass_b[172], pass_cs[6][1], final_R[179]);
	PE2(23, pass_a[117], pass_b[122], pass_cs[5][2], pass_a[143], pass_b[148], pass_cs[5][3], final_R[154]);
	PE2(22, pass_a[92], pass_b[96], pass_cs[4][4], pass_a[119], pass_b[123], pass_cs[4][5], final_R[128]);
	PE2(21, pass_a[66], pass_b[69], pass_cs[3][6], pass_a[94], pass_b[97], pass_cs[3][7], final_R[101]);
	PE2(20, pass_a[39], pass_b[41], pass_cs[2][8], pass_a[68], pass_b[70], pass_cs[2][9], final_R[73]);
	PE2(19, pass_a[11], pass_b[12], pass_cs[1][10], pass_a[41], pass_b[42], pass_cs[1][11], final_R[44]);
	PE2(18, Feed_A[14], Feed_B[14], pass_cs[0][12], pass_a[13], pass_b[13], pass_cs[0][13], final_R[14]);

	PE3(7,out_cs_Q[7],Q_left[6], Q_right[6],Q_left[7], Q_right[7], final_Q_left[7], final_Q_right[7]);
	PE2(24, pass_a[165], pass_b[172], PE1_pass_cs[7], pass_a[189], pass_b[196], pass_cs[7][0], final_R[204]);
	PE2(23, pass_a[142], pass_b[148], pass_cs[6][1], pass_a[167], pass_b[173], pass_cs[6][2], final_R[180]);
	PE2(22, pass_a[118], pass_b[123], pass_cs[5][3], pass_a[144], pass_b[149], pass_cs[5][4], final_R[155]);
	PE2(21, pass_a[93], pass_b[97], pass_cs[4][5], pass_a[120], pass_b[124], pass_cs[4][6], final_R[129]);
	PE2(20, pass_a[67], pass_b[70], pass_cs[3][7], pass_a[95], pass_b[98], pass_cs[3][8], final_R[102]);
	PE2(19, pass_a[40], pass_b[42], pass_cs[2][9], pass_a[69], pass_b[71], pass_cs[2][10], final_R[74]);
	PE2(18, pass_a[12], pass_b[13], pass_cs[1][11], pass_a[42], pass_b[43], pass_cs[1][12], final_R[45]);
	PE2(17, Feed_A[15], Feed_B[15], pass_cs[0][13], pass_a[14], pass_b[14], pass_cs[0][14], final_R[15]);

	PE1(24, out_a[7], pass_b[196], out_a[8], PE1_pass_cs[8], out_cs_Q[8], final_R[228]);
	PE2(23, pass_a[166], pass_b[173], pass_cs[7][0], pass_a[190], pass_b[197], pass_cs[7][1], final_R[205]);
	PE2(22, pass_a[143], pass_b[149], pass_cs[6][2], pass_a[168], pass_b[174], pass_cs[6][3], final_R[181]);
	PE2(21, pass_a[119], pass_b[124], pass_cs[5][4], pass_a[145], pass_b[150], pass_cs[5][5], final_R[156]);
	PE2(20, pass_a[94], pass_b[98], pass_cs[4][6], pass_a[121], pass_b[125], pass_cs[4][7], final_R[130]);
	PE2(19, pass_a[68], pass_b[71], pass_cs[3][8], pass_a[96], pass_b[99], pass_cs[3][9], final_R[103]);
	PE2(18, pass_a[41], pass_b[43], pass_cs[2][10], pass_a[70], pass_b[72], pass_cs[2][11], final_R[75]);
	PE2(17, pass_a[13], pass_b[14], pass_cs[1][12], pass_a[43], pass_b[44], pass_cs[1][13], final_R[46]);
	PE2(16, Feed_A[16], Feed_B[16], pass_cs[0][14], pass_a[15], pass_b[15], pass_cs[0][15], final_R[16]);

	PE3(8,out_cs_Q[8],Q_left[7], Q_right[7],Q_left[8], Q_right[8], final_Q_left[8], final_Q_right[8]);
	PE2(23, pass_a[189], pass_b[197], PE1_pass_cs[8], pass_a[212], pass_b[220], pass_cs[8][0], final_R[229]);
	PE2(22, pass_a[167], pass_b[174], pass_cs[7][1], pass_a[191], pass_b[198], pass_cs[7][2], final_R[206]);
	PE2(21, pass_a[144], pass_b[150], pass_cs[6][3], pass_a[169], pass_b[175], pass_cs[6][4], final_R[182]);
	PE2(20, pass_a[120], pass_b[125], pass_cs[5][5], pass_a[146], pass_b[151], pass_cs[5][6], final_R[157]);
	PE2(19, pass_a[95], pass_b[99], pass_cs[4][7], pass_a[122], pass_b[126], pass_cs[4][8], final_R[131]);
	PE2(18, pass_a[69], pass_b[72], pass_cs[3][9], pass_a[97], pass_b[100], pass_cs[3][10], final_R[104]);
	PE2(17, pass_a[42], pass_b[44], pass_cs[2][11], pass_a[71], pass_b[73], pass_cs[2][12], final_R[76]);
	PE2(16, pass_a[14], pass_b[15], pass_cs[1][13], pass_a[44], pass_b[45], pass_cs[1][14], final_R[47]);
	PE2(15, Feed_A[17], Feed_B[17], pass_cs[0][15], pass_a[16], pass_b[16], pass_cs[0][16], final_R[17]);

	PE1(23, out_a[8], pass_b[220], out_a[9], PE1_pass_cs[9], out_cs_Q[9], final_R[252]);
	PE2(22, pass_a[190], pass_b[198], pass_cs[8][0], pass_a[213], pass_b[221], pass_cs[8][1], final_R[230]);
	PE2(21, pass_a[168], pass_b[175], pass_cs[7][2], pass_a[192], pass_b[199], pass_cs[7][3], final_R[207]);
	PE2(20, pass_a[145], pass_b[151], pass_cs[6][4], pass_a[170], pass_b[176], pass_cs[6][5], final_R[183]);
	PE2(19, pass_a[121], pass_b[126], pass_cs[5][6], pass_a[147], pass_b[152], pass_cs[5][7], final_R[158]);
	PE2(18, pass_a[96], pass_b[100], pass_cs[4][8], pass_a[123], pass_b[127], pass_cs[4][9], final_R[132]);
	PE2(17, pass_a[70], pass_b[73], pass_cs[3][10], pass_a[98], pass_b[101], pass_cs[3][11], final_R[105]);
	PE2(16, pass_a[43], pass_b[45], pass_cs[2][12], pass_a[72], pass_b[74], pass_cs[2][13], final_R[77]);
	PE2(15, pass_a[15], pass_b[16], pass_cs[1][14], pass_a[45], pass_b[46], pass_cs[1][15], final_R[48]);
	PE2(14, Feed_A[18], Feed_B[18], pass_cs[0][16], pass_a[17], pass_b[17], pass_cs[0][17], final_R[18]);

	PE3(9,out_cs_Q[9],Q_left[8], Q_right[8],Q_left[9], Q_right[9], final_Q_left[9], final_Q_right[9]);
	PE2(22, pass_a[212], pass_b[221], PE1_pass_cs[9], pass_a[234], pass_b[243], pass_cs[9][0], final_R[253]);
	PE2(21, pass_a[191], pass_b[199], pass_cs[8][1], pass_a[214], pass_b[222], pass_cs[8][2], final_R[231]);
	PE2(20, pass_a[169], pass_b[176], pass_cs[7][3], pass_a[193], pass_b[200], pass_cs[7][4], final_R[208]);
	PE2(19, pass_a[146], pass_b[152], pass_cs[6][5], pass_a[171], pass_b[177], pass_cs[6][6], final_R[184]);
	PE2(18, pass_a[122], pass_b[127], pass_cs[5][7], pass_a[148], pass_b[153], pass_cs[5][8], final_R[159]);
	PE2(17, pass_a[97], pass_b[101], pass_cs[4][9], pass_a[124], pass_b[128], pass_cs[4][10], final_R[133]);
	PE2(16, pass_a[71], pass_b[74], pass_cs[3][11], pass_a[99], pass_b[102], pass_cs[3][12], final_R[106]);
	PE2(15, pass_a[44], pass_b[46], pass_cs[2][13], pass_a[73], pass_b[75], pass_cs[2][14], final_R[78]);
	PE2(14, pass_a[16], pass_b[17], pass_cs[1][15], pass_a[46], pass_b[47], pass_cs[1][16], final_R[49]);
	PE2(13, Feed_A[19], Feed_B[19], pass_cs[0][17], pass_a[18], pass_b[18], pass_cs[0][18], final_R[19]);

	PE1(22, out_a[9], pass_b[243], out_a[10], PE1_pass_cs[10], out_cs_Q[10], final_R[275]);
	PE2(21, pass_a[213], pass_b[222], pass_cs[9][0], pass_a[235], pass_b[244], pass_cs[9][1], final_R[254]);
	PE2(20, pass_a[192], pass_b[200], pass_cs[8][2], pass_a[215], pass_b[223], pass_cs[8][3], final_R[232]);
	PE2(19, pass_a[170], pass_b[177], pass_cs[7][4], pass_a[194], pass_b[201], pass_cs[7][5], final_R[209]);
	PE2(18, pass_a[147], pass_b[153], pass_cs[6][6], pass_a[172], pass_b[178], pass_cs[6][7], final_R[185]);
	PE2(17, pass_a[123], pass_b[128], pass_cs[5][8], pass_a[149], pass_b[154], pass_cs[5][9], final_R[160]);
	PE2(16, pass_a[98], pass_b[102], pass_cs[4][10], pass_a[125], pass_b[129], pass_cs[4][11], final_R[134]);
	PE2(15, pass_a[72], pass_b[75], pass_cs[3][12], pass_a[100], pass_b[103], pass_cs[3][13], final_R[107]);
	PE2(14, pass_a[45], pass_b[47], pass_cs[2][14], pass_a[74], pass_b[76], pass_cs[2][15], final_R[79]);
	PE2(13, pass_a[17], pass_b[18], pass_cs[1][16], pass_a[47], pass_b[48], pass_cs[1][17], final_R[50]);
	PE2(12, Feed_A[20], Feed_B[20], pass_cs[0][18], pass_a[19], pass_b[19], pass_cs[0][19], final_R[20]);

	PE3(10,out_cs_Q[10],Q_left[9], Q_right[9],Q_left[10], Q_right[10], final_Q_left[10], final_Q_right[10]);
	PE2(21, pass_a[234], pass_b[244], PE1_pass_cs[10], pass_a[255], pass_b[265], pass_cs[10][0], final_R[276]);
	PE2(20, pass_a[214], pass_b[223], pass_cs[9][1], pass_a[236], pass_b[245], pass_cs[9][2], final_R[255]);
	PE2(19, pass_a[193], pass_b[201], pass_cs[8][3], pass_a[216], pass_b[224], pass_cs[8][4], final_R[233]);
	PE2(18, pass_a[171], pass_b[178], pass_cs[7][5], pass_a[195], pass_b[202], pass_cs[7][6], final_R[210]);
	PE2(17, pass_a[148], pass_b[154], pass_cs[6][7], pass_a[173], pass_b[179], pass_cs[6][8], final_R[186]);
	PE2(16, pass_a[124], pass_b[129], pass_cs[5][9], pass_a[150], pass_b[155], pass_cs[5][10], final_R[161]);
	PE2(15, pass_a[99], pass_b[103], pass_cs[4][11], pass_a[126], pass_b[130], pass_cs[4][12], final_R[135]);
	PE2(14, pass_a[73], pass_b[76], pass_cs[3][13], pass_a[101], pass_b[104], pass_cs[3][14], final_R[108]);
	PE2(13, pass_a[46], pass_b[48], pass_cs[2][15], pass_a[75], pass_b[77], pass_cs[2][16], final_R[80]);
	PE2(12, pass_a[18], pass_b[19], pass_cs[1][17], pass_a[48], pass_b[49], pass_cs[1][18], final_R[51]);
	PE2(11, Feed_A[21], Feed_B[21], pass_cs[0][19], pass_a[20], pass_b[20], pass_cs[0][20], final_R[21]);

	PE1(21, out_a[10], pass_b[265], out_a[11], PE1_pass_cs[11], out_cs_Q[11], final_R[297]);
	PE2(20, pass_a[235], pass_b[245], pass_cs[10][0], pass_a[256], pass_b[266], pass_cs[10][1], final_R[277]);
	PE2(19, pass_a[215], pass_b[224], pass_cs[9][2], pass_a[237], pass_b[246], pass_cs[9][3], final_R[256]);
	PE2(18, pass_a[194], pass_b[202], pass_cs[8][4], pass_a[217], pass_b[225], pass_cs[8][5], final_R[234]);
	PE2(17, pass_a[172], pass_b[179], pass_cs[7][6], pass_a[196], pass_b[203], pass_cs[7][7], final_R[211]);
	PE2(16, pass_a[149], pass_b[155], pass_cs[6][8], pass_a[174], pass_b[180], pass_cs[6][9], final_R[187]);
	PE2(15, pass_a[125], pass_b[130], pass_cs[5][10], pass_a[151], pass_b[156], pass_cs[5][11], final_R[162]);
	PE2(14, pass_a[100], pass_b[104], pass_cs[4][12], pass_a[127], pass_b[131], pass_cs[4][13], final_R[136]);
	PE2(13, pass_a[74], pass_b[77], pass_cs[3][14], pass_a[102], pass_b[105], pass_cs[3][15], final_R[109]);
	PE2(12, pass_a[47], pass_b[49], pass_cs[2][16], pass_a[76], pass_b[78], pass_cs[2][17], final_R[81]);
	PE2(11, pass_a[19], pass_b[20], pass_cs[1][18], pass_a[49], pass_b[50], pass_cs[1][19], final_R[52]);
	PE2(10, Feed_A[22], Feed_B[22], pass_cs[0][20], pass_a[21], pass_b[21], pass_cs[0][21], final_R[22]);

	PE3(11,out_cs_Q[11],Q_left[10], Q_right[10],Q_left[11], Q_right[11], final_Q_left[11], final_Q_right[11]);
	PE2(20, pass_a[255], pass_b[266], PE1_pass_cs[11], pass_a[275], pass_b[286], pass_cs[11][0], final_R[298]);
	PE2(19, pass_a[236], pass_b[246], pass_cs[10][1], pass_a[257], pass_b[267], pass_cs[10][2], final_R[278]);
	PE2(18, pass_a[216], pass_b[225], pass_cs[9][3], pass_a[238], pass_b[247], pass_cs[9][4], final_R[257]);
	PE2(17, pass_a[195], pass_b[203], pass_cs[8][5], pass_a[218], pass_b[226], pass_cs[8][6], final_R[235]);
	PE2(16, pass_a[173], pass_b[180], pass_cs[7][7], pass_a[197], pass_b[204], pass_cs[7][8], final_R[212]);
	PE2(15, pass_a[150], pass_b[156], pass_cs[6][9], pass_a[175], pass_b[181], pass_cs[6][10], final_R[188]);
	PE2(14, pass_a[126], pass_b[131], pass_cs[5][11], pass_a[152], pass_b[157], pass_cs[5][12], final_R[163]);
	PE2(13, pass_a[101], pass_b[105], pass_cs[4][13], pass_a[128], pass_b[132], pass_cs[4][14], final_R[137]);
	PE2(12, pass_a[75], pass_b[78], pass_cs[3][15], pass_a[103], pass_b[106], pass_cs[3][16], final_R[110]);
	PE2(11, pass_a[48], pass_b[50], pass_cs[2][17], pass_a[77], pass_b[79], pass_cs[2][18], final_R[82]);
	PE2(10, pass_a[20], pass_b[21], pass_cs[1][19], pass_a[50], pass_b[51], pass_cs[1][20], final_R[53]);
	PE2(9, Feed_A[23], Feed_B[23], pass_cs[0][21], pass_a[22], pass_b[22], pass_cs[0][22], final_R[23]);

	PE1(20, out_a[11], pass_b[286], out_a[12], PE1_pass_cs[12], out_cs_Q[12], final_R[318]);
	PE2(19, pass_a[256], pass_b[267], pass_cs[11][0], pass_a[276], pass_b[287], pass_cs[11][1], final_R[299]);
	PE2(18, pass_a[237], pass_b[247], pass_cs[10][2], pass_a[258], pass_b[268], pass_cs[10][3], final_R[279]);
	PE2(17, pass_a[217], pass_b[226], pass_cs[9][4], pass_a[239], pass_b[248], pass_cs[9][5], final_R[258]);
	PE2(16, pass_a[196], pass_b[204], pass_cs[8][6], pass_a[219], pass_b[227], pass_cs[8][7], final_R[236]);
	PE2(15, pass_a[174], pass_b[181], pass_cs[7][8], pass_a[198], pass_b[205], pass_cs[7][9], final_R[213]);
	PE2(14, pass_a[151], pass_b[157], pass_cs[6][10], pass_a[176], pass_b[182], pass_cs[6][11], final_R[189]);
	PE2(13, pass_a[127], pass_b[132], pass_cs[5][12], pass_a[153], pass_b[158], pass_cs[5][13], final_R[164]);
	PE2(12, pass_a[102], pass_b[106], pass_cs[4][14], pass_a[129], pass_b[133], pass_cs[4][15], final_R[138]);
	PE2(11, pass_a[76], pass_b[79], pass_cs[3][16], pass_a[104], pass_b[107], pass_cs[3][17], final_R[111]);
	PE2(10, pass_a[49], pass_b[51], pass_cs[2][18], pass_a[78], pass_b[80], pass_cs[2][19], final_R[83]);
	PE2(9, pass_a[21], pass_b[22], pass_cs[1][20], pass_a[51], pass_b[52], pass_cs[1][21], final_R[54]);
	PE2(8, Feed_A[24], Feed_B[24], pass_cs[0][22], pass_a[23], pass_b[23], pass_cs[0][23], final_R[24]);

	PE3(12,out_cs_Q[12],Q_left[11], Q_right[11],Q_left[12], Q_right[12], final_Q_left[12], final_Q_right[12]);
	PE2(19, pass_a[275], pass_b[287], PE1_pass_cs[12], pass_a[294], pass_b[306], pass_cs[12][0], final_R[319]);
	PE2(18, pass_a[257], pass_b[268], pass_cs[11][1], pass_a[277], pass_b[288], pass_cs[11][2], final_R[300]);
	PE2(17, pass_a[238], pass_b[248], pass_cs[10][3], pass_a[259], pass_b[269], pass_cs[10][4], final_R[280]);
	PE2(16, pass_a[218], pass_b[227], pass_cs[9][5], pass_a[240], pass_b[249], pass_cs[9][6], final_R[259]);
	PE2(15, pass_a[197], pass_b[205], pass_cs[8][7], pass_a[220], pass_b[228], pass_cs[8][8], final_R[237]);
	PE2(14, pass_a[175], pass_b[182], pass_cs[7][9], pass_a[199], pass_b[206], pass_cs[7][10], final_R[214]);
	PE2(13, pass_a[152], pass_b[158], pass_cs[6][11], pass_a[177], pass_b[183], pass_cs[6][12], final_R[190]);
	PE2(12, pass_a[128], pass_b[133], pass_cs[5][13], pass_a[154], pass_b[159], pass_cs[5][14], final_R[165]);
	PE2(11, pass_a[103], pass_b[107], pass_cs[4][15], pass_a[130], pass_b[134], pass_cs[4][16], final_R[139]);
	PE2(10, pass_a[77], pass_b[80], pass_cs[3][17], pass_a[105], pass_b[108], pass_cs[3][18], final_R[112]);
	PE2(9, pass_a[50], pass_b[52], pass_cs[2][19], pass_a[79], pass_b[81], pass_cs[2][20], final_R[84]);
	PE2(8, pass_a[22], pass_b[23], pass_cs[1][21], pass_a[52], pass_b[53], pass_cs[1][22], final_R[55]);
	PE2(7, Feed_A[25], Feed_B[25], pass_cs[0][23], pass_a[24], pass_b[24], pass_cs[0][24], final_R[25]);

	PE1(19, out_a[12], pass_b[306], out_a[13], PE1_pass_cs[13], out_cs_Q[13], final_R[338]);
	PE2(18, pass_a[276], pass_b[288], pass_cs[12][0], pass_a[295], pass_b[307], pass_cs[12][1], final_R[320]);
	PE2(17, pass_a[258], pass_b[269], pass_cs[11][2], pass_a[278], pass_b[289], pass_cs[11][3], final_R[301]);
	PE2(16, pass_a[239], pass_b[249], pass_cs[10][4], pass_a[260], pass_b[270], pass_cs[10][5], final_R[281]);
	PE2(15, pass_a[219], pass_b[228], pass_cs[9][6], pass_a[241], pass_b[250], pass_cs[9][7], final_R[260]);
	PE2(14, pass_a[198], pass_b[206], pass_cs[8][8], pass_a[221], pass_b[229], pass_cs[8][9], final_R[238]);
	PE2(13, pass_a[176], pass_b[183], pass_cs[7][10], pass_a[200], pass_b[207], pass_cs[7][11], final_R[215]);
	PE2(12, pass_a[153], pass_b[159], pass_cs[6][12], pass_a[178], pass_b[184], pass_cs[6][13], final_R[191]);
	PE2(11, pass_a[129], pass_b[134], pass_cs[5][14], pass_a[155], pass_b[160], pass_cs[5][15], final_R[166]);
	PE2(10, pass_a[104], pass_b[108], pass_cs[4][16], pass_a[131], pass_b[135], pass_cs[4][17], final_R[140]);
	PE2(9, pass_a[78], pass_b[81], pass_cs[3][18], pass_a[106], pass_b[109], pass_cs[3][19], final_R[113]);
	PE2(8, pass_a[51], pass_b[53], pass_cs[2][20], pass_a[80], pass_b[82], pass_cs[2][21], final_R[85]);
	PE2(7, pass_a[23], pass_b[24], pass_cs[1][22], pass_a[53], pass_b[54], pass_cs[1][23], final_R[56]);
	PE2(6, Feed_A[26], Feed_B[26], pass_cs[0][24], pass_a[25], pass_b[25], pass_cs[0][25], final_R[26]);

	PE3(13,out_cs_Q[13],Q_left[12], Q_right[12],Q_left[13], Q_right[13], final_Q_left[13], final_Q_right[13]);
	PE2(18, pass_a[294], pass_b[307], PE1_pass_cs[13], pass_a[312], pass_b[325], pass_cs[13][0], final_R[339]);
	PE2(17, pass_a[277], pass_b[289], pass_cs[12][1], pass_a[296], pass_b[308], pass_cs[12][2], final_R[321]);
	PE2(16, pass_a[259], pass_b[270], pass_cs[11][3], pass_a[279], pass_b[290], pass_cs[11][4], final_R[302]);
	PE2(15, pass_a[240], pass_b[250], pass_cs[10][5], pass_a[261], pass_b[271], pass_cs[10][6], final_R[282]);
	PE2(14, pass_a[220], pass_b[229], pass_cs[9][7], pass_a[242], pass_b[251], pass_cs[9][8], final_R[261]);
	PE2(13, pass_a[199], pass_b[207], pass_cs[8][9], pass_a[222], pass_b[230], pass_cs[8][10], final_R[239]);
	PE2(12, pass_a[177], pass_b[184], pass_cs[7][11], pass_a[201], pass_b[208], pass_cs[7][12], final_R[216]);
	PE2(11, pass_a[154], pass_b[160], pass_cs[6][13], pass_a[179], pass_b[185], pass_cs[6][14], final_R[192]);
	PE2(10, pass_a[130], pass_b[135], pass_cs[5][15], pass_a[156], pass_b[161], pass_cs[5][16], final_R[167]);
	PE2(9, pass_a[105], pass_b[109], pass_cs[4][17], pass_a[132], pass_b[136], pass_cs[4][18], final_R[141]);
	PE2(8, pass_a[79], pass_b[82], pass_cs[3][19], pass_a[107], pass_b[110], pass_cs[3][20], final_R[114]);
	PE2(7, pass_a[52], pass_b[54], pass_cs[2][21], pass_a[81], pass_b[83], pass_cs[2][22], final_R[86]);
	PE2(6, pass_a[24], pass_b[25], pass_cs[1][23], pass_a[54], pass_b[55], pass_cs[1][24], final_R[57]);
	PE2(5, Feed_A[27], Feed_B[27], pass_cs[0][25], pass_a[26], pass_b[26], pass_cs[0][26], final_R[27]);

	PE1(18, out_a[13], pass_b[325], out_a[14], PE1_pass_cs[14], out_cs_Q[14], final_R[357]);
	PE2(17, pass_a[295], pass_b[308], pass_cs[13][0], pass_a[313], pass_b[326], pass_cs[13][1], final_R[340]);
	PE2(16, pass_a[278], pass_b[290], pass_cs[12][2], pass_a[297], pass_b[309], pass_cs[12][3], final_R[322]);
	PE2(15, pass_a[260], pass_b[271], pass_cs[11][4], pass_a[280], pass_b[291], pass_cs[11][5], final_R[303]);
	PE2(14, pass_a[241], pass_b[251], pass_cs[10][6], pass_a[262], pass_b[272], pass_cs[10][7], final_R[283]);
	PE2(13, pass_a[221], pass_b[230], pass_cs[9][8], pass_a[243], pass_b[252], pass_cs[9][9], final_R[262]);
	PE2(12, pass_a[200], pass_b[208], pass_cs[8][10], pass_a[223], pass_b[231], pass_cs[8][11], final_R[240]);
	PE2(11, pass_a[178], pass_b[185], pass_cs[7][12], pass_a[202], pass_b[209], pass_cs[7][13], final_R[217]);
	PE2(10, pass_a[155], pass_b[161], pass_cs[6][14], pass_a[180], pass_b[186], pass_cs[6][15], final_R[193]);
	PE2(9, pass_a[131], pass_b[136], pass_cs[5][16], pass_a[157], pass_b[162], pass_cs[5][17], final_R[168]);
	PE2(8, pass_a[106], pass_b[110], pass_cs[4][18], pass_a[133], pass_b[137], pass_cs[4][19], final_R[142]);
	PE2(7, pass_a[80], pass_b[83], pass_cs[3][20], pass_a[108], pass_b[111], pass_cs[3][21], final_R[115]);
	PE2(6, pass_a[53], pass_b[55], pass_cs[2][22], pass_a[82], pass_b[84], pass_cs[2][23], final_R[87]);
	PE2(5, pass_a[25], pass_b[26], pass_cs[1][24], pass_a[55], pass_b[56], pass_cs[1][25], final_R[58]);
	PE2(4, Feed_A[28], Feed_B[28], pass_cs[0][26], pass_a[27], pass_b[27], pass_cs[0][27], final_R[28]);

	PE3(14,out_cs_Q[14],Q_left[13], Q_right[13],Q_left[14], Q_right[14], final_Q_left[14], final_Q_right[14]);
	PE2(17, pass_a[312], pass_b[326], PE1_pass_cs[14], pass_a[329], pass_b[343], pass_cs[14][0], final_R[358]);
	PE2(16, pass_a[296], pass_b[309], pass_cs[13][1], pass_a[314], pass_b[327], pass_cs[13][2], final_R[341]);
	PE2(15, pass_a[279], pass_b[291], pass_cs[12][3], pass_a[298], pass_b[310], pass_cs[12][4], final_R[323]);
	PE2(14, pass_a[261], pass_b[272], pass_cs[11][5], pass_a[281], pass_b[292], pass_cs[11][6], final_R[304]);
	PE2(13, pass_a[242], pass_b[252], pass_cs[10][7], pass_a[263], pass_b[273], pass_cs[10][8], final_R[284]);
	PE2(12, pass_a[222], pass_b[231], pass_cs[9][9], pass_a[244], pass_b[253], pass_cs[9][10], final_R[263]);
	PE2(11, pass_a[201], pass_b[209], pass_cs[8][11], pass_a[224], pass_b[232], pass_cs[8][12], final_R[241]);
	PE2(10, pass_a[179], pass_b[186], pass_cs[7][13], pass_a[203], pass_b[210], pass_cs[7][14], final_R[218]);
	PE2(9, pass_a[156], pass_b[162], pass_cs[6][15], pass_a[181], pass_b[187], pass_cs[6][16], final_R[194]);
	PE2(8, pass_a[132], pass_b[137], pass_cs[5][17], pass_a[158], pass_b[163], pass_cs[5][18], final_R[169]);
	PE2(7, pass_a[107], pass_b[111], pass_cs[4][19], pass_a[134], pass_b[138], pass_cs[4][20], final_R[143]);
	PE2(6, pass_a[81], pass_b[84], pass_cs[3][21], pass_a[109], pass_b[112], pass_cs[3][22], final_R[116]);
	PE2(5, pass_a[54], pass_b[56], pass_cs[2][23], pass_a[83], pass_b[85], pass_cs[2][24], final_R[88]);
	PE2(4, pass_a[26], pass_b[27], pass_cs[1][25], pass_a[56], pass_b[57], pass_cs[1][26], final_R[59]);
	PE2(3, Feed_A[29], Feed_B[29], pass_cs[0][27], pass_a[28], pass_b[28], pass_cs[0][28], final_R[29]);

	PE1(17, out_a[14], pass_b[343], out_a[15], PE1_pass_cs[15], out_cs_Q[15], final_R[375]);
	PE2(16, pass_a[313], pass_b[327], pass_cs[14][0], pass_a[330], pass_b[344], pass_cs[14][1], final_R[359]);
	PE2(15, pass_a[297], pass_b[310], pass_cs[13][2], pass_a[315], pass_b[328], pass_cs[13][3], final_R[342]);
	PE2(14, pass_a[280], pass_b[292], pass_cs[12][4], pass_a[299], pass_b[311], pass_cs[12][5], final_R[324]);
	PE2(13, pass_a[262], pass_b[273], pass_cs[11][6], pass_a[282], pass_b[293], pass_cs[11][7], final_R[305]);
	PE2(12, pass_a[243], pass_b[253], pass_cs[10][8], pass_a[264], pass_b[274], pass_cs[10][9], final_R[285]);
	PE2(11, pass_a[223], pass_b[232], pass_cs[9][10], pass_a[245], pass_b[254], pass_cs[9][11], final_R[264]);
	PE2(10, pass_a[202], pass_b[210], pass_cs[8][12], pass_a[225], pass_b[233], pass_cs[8][13], final_R[242]);
	PE2(9, pass_a[180], pass_b[187], pass_cs[7][14], pass_a[204], pass_b[211], pass_cs[7][15], final_R[219]);
	PE2(8, pass_a[157], pass_b[163], pass_cs[6][16], pass_a[182], pass_b[188], pass_cs[6][17], final_R[195]);
	PE2(7, pass_a[133], pass_b[138], pass_cs[5][18], pass_a[159], pass_b[164], pass_cs[5][19], final_R[170]);
	PE2(6, pass_a[108], pass_b[112], pass_cs[4][20], pass_a[135], pass_b[139], pass_cs[4][21], final_R[144]);
	PE2(5, pass_a[82], pass_b[85], pass_cs[3][22], pass_a[110], pass_b[113], pass_cs[3][23], final_R[117]);
	PE2(4, pass_a[55], pass_b[57], pass_cs[2][24], pass_a[84], pass_b[86], pass_cs[2][25], final_R[89]);
	PE2(3, pass_a[27], pass_b[28], pass_cs[1][26], pass_a[57], pass_b[58], pass_cs[1][27], final_R[60]);
	PE2(2, Feed_A[30], Feed_B[30], pass_cs[0][28], pass_a[29], pass_b[29], pass_cs[0][29], final_R[30]);

	PE3(15,out_cs_Q[15],Q_left[14], Q_right[14],Q_left[15], Q_right[15], final_Q_left[15], final_Q_right[15]);
	PE2(16, pass_a[329], pass_b[344], PE1_pass_cs[15], pass_a[345], pass_b[360], pass_cs[15][0], final_R[376]);
	PE2(15, pass_a[314], pass_b[328], pass_cs[14][1], pass_a[331], pass_b[345], pass_cs[14][2], final_R[360]);
	PE2(14, pass_a[298], pass_b[311], pass_cs[13][3], pass_a[316], pass_b[329], pass_cs[13][4], final_R[343]);
	PE2(13, pass_a[281], pass_b[293], pass_cs[12][5], pass_a[300], pass_b[312], pass_cs[12][6], final_R[325]);
	PE2(12, pass_a[263], pass_b[274], pass_cs[11][7], pass_a[283], pass_b[294], pass_cs[11][8], final_R[306]);
	PE2(11, pass_a[244], pass_b[254], pass_cs[10][9], pass_a[265], pass_b[275], pass_cs[10][10], final_R[286]);
	PE2(10, pass_a[224], pass_b[233], pass_cs[9][11], pass_a[246], pass_b[255], pass_cs[9][12], final_R[265]);
	PE2(9, pass_a[203], pass_b[211], pass_cs[8][13], pass_a[226], pass_b[234], pass_cs[8][14], final_R[243]);
	PE2(8, pass_a[181], pass_b[188], pass_cs[7][15], pass_a[205], pass_b[212], pass_cs[7][16], final_R[220]);
	PE2(7, pass_a[158], pass_b[164], pass_cs[6][17], pass_a[183], pass_b[189], pass_cs[6][18], final_R[196]);
	PE2(6, pass_a[134], pass_b[139], pass_cs[5][19], pass_a[160], pass_b[165], pass_cs[5][20], final_R[171]);
	PE2(5, pass_a[109], pass_b[113], pass_cs[4][21], pass_a[136], pass_b[140], pass_cs[4][22], final_R[145]);
	PE2(4, pass_a[83], pass_b[86], pass_cs[3][23], pass_a[111], pass_b[114], pass_cs[3][24], final_R[118]);
	PE2(3, pass_a[56], pass_b[58], pass_cs[2][25], pass_a[85], pass_b[87], pass_cs[2][26], final_R[90]);
	PE2(2, pass_a[28], pass_b[29], pass_cs[1][27], pass_a[58], pass_b[59], pass_cs[1][28], final_R[61]);
	PE2_tail(Feed_A[31], Feed_B[31], pass_cs[0][29], pass_b[30], final_R[31]);

	PE1(16, out_a[15], pass_b[360], out_a[16], PE1_pass_cs[16], out_cs_Q[16], final_R[392]);
	PE2(15, pass_a[330], pass_b[345], pass_cs[15][0], pass_a[346], pass_b[361], pass_cs[15][1], final_R[377]);
	PE2(14, pass_a[315], pass_b[329], pass_cs[14][2], pass_a[332], pass_b[346], pass_cs[14][3], final_R[361]);
	PE2(13, pass_a[299], pass_b[312], pass_cs[13][4], pass_a[317], pass_b[330], pass_cs[13][5], final_R[344]);
	PE2(12, pass_a[282], pass_b[294], pass_cs[12][6], pass_a[301], pass_b[313], pass_cs[12][7], final_R[326]);
	PE2(11, pass_a[264], pass_b[275], pass_cs[11][8], pass_a[284], pass_b[295], pass_cs[11][9], final_R[307]);
	PE2(10, pass_a[245], pass_b[255], pass_cs[10][10], pass_a[266], pass_b[276], pass_cs[10][11], final_R[287]);
	PE2(9, pass_a[225], pass_b[234], pass_cs[9][12], pass_a[247], pass_b[256], pass_cs[9][13], final_R[266]);
	PE2(8, pass_a[204], pass_b[212], pass_cs[8][14], pass_a[227], pass_b[235], pass_cs[8][15], final_R[244]);
	PE2(7, pass_a[182], pass_b[189], pass_cs[7][16], pass_a[206], pass_b[213], pass_cs[7][17], final_R[221]);
	PE2(6, pass_a[159], pass_b[165], pass_cs[6][18], pass_a[184], pass_b[190], pass_cs[6][19], final_R[197]);
	PE2(5, pass_a[135], pass_b[140], pass_cs[5][20], pass_a[161], pass_b[166], pass_cs[5][21], final_R[172]);
	PE2(4, pass_a[110], pass_b[114], pass_cs[4][22], pass_a[137], pass_b[141], pass_cs[4][23], final_R[146]);
	PE2(3, pass_a[84], pass_b[87], pass_cs[3][24], pass_a[112], pass_b[115], pass_cs[3][25], final_R[119]);
	PE2(2, pass_a[57], pass_b[59], pass_cs[2][26], pass_a[86], pass_b[88], pass_cs[2][27], final_R[91]);
	PE2_tail(pass_a[29], pass_b[30], pass_cs[1][28], pass_b[60], final_R[62]);

	PE3(16,out_cs_Q[16],Q_left[15], Q_right[15],Q_left[16], Q_right[16], final_Q_left[16], final_Q_right[16]);
	PE2(15, pass_a[345], pass_b[361], PE1_pass_cs[16], pass_a[360], pass_b[376], pass_cs[16][0], final_R[393]);
	PE2(14, pass_a[331], pass_b[346], pass_cs[15][1], pass_a[347], pass_b[362], pass_cs[15][2], final_R[378]);
	PE2(13, pass_a[316], pass_b[330], pass_cs[14][3], pass_a[333], pass_b[347], pass_cs[14][4], final_R[362]);
	PE2(12, pass_a[300], pass_b[313], pass_cs[13][5], pass_a[318], pass_b[331], pass_cs[13][6], final_R[345]);
	PE2(11, pass_a[283], pass_b[295], pass_cs[12][7], pass_a[302], pass_b[314], pass_cs[12][8], final_R[327]);
	PE2(10, pass_a[265], pass_b[276], pass_cs[11][9], pass_a[285], pass_b[296], pass_cs[11][10], final_R[308]);
	PE2(9, pass_a[246], pass_b[256], pass_cs[10][11], pass_a[267], pass_b[277], pass_cs[10][12], final_R[288]);
	PE2(8, pass_a[226], pass_b[235], pass_cs[9][13], pass_a[248], pass_b[257], pass_cs[9][14], final_R[267]);
	PE2(7, pass_a[205], pass_b[213], pass_cs[8][15], pass_a[228], pass_b[236], pass_cs[8][16], final_R[245]);
	PE2(6, pass_a[183], pass_b[190], pass_cs[7][17], pass_a[207], pass_b[214], pass_cs[7][18], final_R[222]);
	PE2(5, pass_a[160], pass_b[166], pass_cs[6][19], pass_a[185], pass_b[191], pass_cs[6][20], final_R[198]);
	PE2(4, pass_a[136], pass_b[141], pass_cs[5][21], pass_a[162], pass_b[167], pass_cs[5][22], final_R[173]);
	PE2(3, pass_a[111], pass_b[115], pass_cs[4][23], pass_a[138], pass_b[142], pass_cs[4][24], final_R[147]);
	PE2(2, pass_a[85], pass_b[88], pass_cs[3][25], pass_a[113], pass_b[116], pass_cs[3][26], final_R[120]);
	PE2_tail(pass_a[58], pass_b[60], pass_cs[2][27], pass_b[89], final_R[92]);

	PE1(15, out_a[16], pass_b[376], out_a[17], PE1_pass_cs[17], out_cs_Q[17], final_R[408]);
	PE2(14, pass_a[346], pass_b[362], pass_cs[16][0], pass_a[361], pass_b[377], pass_cs[16][1], final_R[394]);
	PE2(13, pass_a[332], pass_b[347], pass_cs[15][2], pass_a[348], pass_b[363], pass_cs[15][3], final_R[379]);
	PE2(12, pass_a[317], pass_b[331], pass_cs[14][4], pass_a[334], pass_b[348], pass_cs[14][5], final_R[363]);
	PE2(11, pass_a[301], pass_b[314], pass_cs[13][6], pass_a[319], pass_b[332], pass_cs[13][7], final_R[346]);
	PE2(10, pass_a[284], pass_b[296], pass_cs[12][8], pass_a[303], pass_b[315], pass_cs[12][9], final_R[328]);
	PE2(9, pass_a[266], pass_b[277], pass_cs[11][10], pass_a[286], pass_b[297], pass_cs[11][11], final_R[309]);
	PE2(8, pass_a[247], pass_b[257], pass_cs[10][12], pass_a[268], pass_b[278], pass_cs[10][13], final_R[289]);
	PE2(7, pass_a[227], pass_b[236], pass_cs[9][14], pass_a[249], pass_b[258], pass_cs[9][15], final_R[268]);
	PE2(6, pass_a[206], pass_b[214], pass_cs[8][16], pass_a[229], pass_b[237], pass_cs[8][17], final_R[246]);
	PE2(5, pass_a[184], pass_b[191], pass_cs[7][18], pass_a[208], pass_b[215], pass_cs[7][19], final_R[223]);
	PE2(4, pass_a[161], pass_b[167], pass_cs[6][20], pass_a[186], pass_b[192], pass_cs[6][21], final_R[199]);
	PE2(3, pass_a[137], pass_b[142], pass_cs[5][22], pass_a[163], pass_b[168], pass_cs[5][23], final_R[174]);
	PE2(2, pass_a[112], pass_b[116], pass_cs[4][24], pass_a[139], pass_b[143], pass_cs[4][25], final_R[148]);
	PE2_tail(pass_a[86], pass_b[89], pass_cs[3][26], pass_b[117], final_R[121]);

	PE3(17,out_cs_Q[17],Q_left[16], Q_right[16],Q_left[17], Q_right[17], final_Q_left[17], final_Q_right[17]);
	PE2(14, pass_a[360], pass_b[377], PE1_pass_cs[17], pass_a[374], pass_b[391], pass_cs[17][0], final_R[409]);
	PE2(13, pass_a[347], pass_b[363], pass_cs[16][1], pass_a[362], pass_b[378], pass_cs[16][2], final_R[395]);
	PE2(12, pass_a[333], pass_b[348], pass_cs[15][3], pass_a[349], pass_b[364], pass_cs[15][4], final_R[380]);
	PE2(11, pass_a[318], pass_b[332], pass_cs[14][5], pass_a[335], pass_b[349], pass_cs[14][6], final_R[364]);
	PE2(10, pass_a[302], pass_b[315], pass_cs[13][7], pass_a[320], pass_b[333], pass_cs[13][8], final_R[347]);
	PE2(9, pass_a[285], pass_b[297], pass_cs[12][9], pass_a[304], pass_b[316], pass_cs[12][10], final_R[329]);
	PE2(8, pass_a[267], pass_b[278], pass_cs[11][11], pass_a[287], pass_b[298], pass_cs[11][12], final_R[310]);
	PE2(7, pass_a[248], pass_b[258], pass_cs[10][13], pass_a[269], pass_b[279], pass_cs[10][14], final_R[290]);
	PE2(6, pass_a[228], pass_b[237], pass_cs[9][15], pass_a[250], pass_b[259], pass_cs[9][16], final_R[269]);
	PE2(5, pass_a[207], pass_b[215], pass_cs[8][17], pass_a[230], pass_b[238], pass_cs[8][18], final_R[247]);
	PE2(4, pass_a[185], pass_b[192], pass_cs[7][19], pass_a[209], pass_b[216], pass_cs[7][20], final_R[224]);
	PE2(3, pass_a[162], pass_b[168], pass_cs[6][21], pass_a[187], pass_b[193], pass_cs[6][22], final_R[200]);
	PE2(2, pass_a[138], pass_b[143], pass_cs[5][23], pass_a[164], pass_b[169], pass_cs[5][24], final_R[175]);
	PE2_tail(pass_a[113], pass_b[117], pass_cs[4][25], pass_b[144], final_R[149]);

	PE1(14, out_a[17], pass_b[391], out_a[18], PE1_pass_cs[18], out_cs_Q[18], final_R[423]);
	PE2(13, pass_a[361], pass_b[378], pass_cs[17][0], pass_a[375], pass_b[392], pass_cs[17][1], final_R[410]);
	PE2(12, pass_a[348], pass_b[364], pass_cs[16][2], pass_a[363], pass_b[379], pass_cs[16][3], final_R[396]);
	PE2(11, pass_a[334], pass_b[349], pass_cs[15][4], pass_a[350], pass_b[365], pass_cs[15][5], final_R[381]);
	PE2(10, pass_a[319], pass_b[333], pass_cs[14][6], pass_a[336], pass_b[350], pass_cs[14][7], final_R[365]);
	PE2(9, pass_a[303], pass_b[316], pass_cs[13][8], pass_a[321], pass_b[334], pass_cs[13][9], final_R[348]);
	PE2(8, pass_a[286], pass_b[298], pass_cs[12][10], pass_a[305], pass_b[317], pass_cs[12][11], final_R[330]);
	PE2(7, pass_a[268], pass_b[279], pass_cs[11][12], pass_a[288], pass_b[299], pass_cs[11][13], final_R[311]);
	PE2(6, pass_a[249], pass_b[259], pass_cs[10][14], pass_a[270], pass_b[280], pass_cs[10][15], final_R[291]);
	PE2(5, pass_a[229], pass_b[238], pass_cs[9][16], pass_a[251], pass_b[260], pass_cs[9][17], final_R[270]);
	PE2(4, pass_a[208], pass_b[216], pass_cs[8][18], pass_a[231], pass_b[239], pass_cs[8][19], final_R[248]);
	PE2(3, pass_a[186], pass_b[193], pass_cs[7][20], pass_a[210], pass_b[217], pass_cs[7][21], final_R[225]);
	PE2(2, pass_a[163], pass_b[169], pass_cs[6][22], pass_a[188], pass_b[194], pass_cs[6][23], final_R[201]);
	PE2_tail(pass_a[139], pass_b[144], pass_cs[5][24], pass_b[170], final_R[176]);

	PE3(18,out_cs_Q[18],Q_left[17], Q_right[17],Q_left[18], Q_right[18], final_Q_left[18], final_Q_right[18]);
	PE2(13, pass_a[374], pass_b[392], PE1_pass_cs[18], pass_a[387], pass_b[405], pass_cs[18][0], final_R[424]);
	PE2(12, pass_a[362], pass_b[379], pass_cs[17][1], pass_a[376], pass_b[393], pass_cs[17][2], final_R[411]);
	PE2(11, pass_a[349], pass_b[365], pass_cs[16][3], pass_a[364], pass_b[380], pass_cs[16][4], final_R[397]);
	PE2(10, pass_a[335], pass_b[350], pass_cs[15][5], pass_a[351], pass_b[366], pass_cs[15][6], final_R[382]);
	PE2(9, pass_a[320], pass_b[334], pass_cs[14][7], pass_a[337], pass_b[351], pass_cs[14][8], final_R[366]);
	PE2(8, pass_a[304], pass_b[317], pass_cs[13][9], pass_a[322], pass_b[335], pass_cs[13][10], final_R[349]);
	PE2(7, pass_a[287], pass_b[299], pass_cs[12][11], pass_a[306], pass_b[318], pass_cs[12][12], final_R[331]);
	PE2(6, pass_a[269], pass_b[280], pass_cs[11][13], pass_a[289], pass_b[300], pass_cs[11][14], final_R[312]);
	PE2(5, pass_a[250], pass_b[260], pass_cs[10][15], pass_a[271], pass_b[281], pass_cs[10][16], final_R[292]);
	PE2(4, pass_a[230], pass_b[239], pass_cs[9][17], pass_a[252], pass_b[261], pass_cs[9][18], final_R[271]);
	PE2(3, pass_a[209], pass_b[217], pass_cs[8][19], pass_a[232], pass_b[240], pass_cs[8][20], final_R[249]);
	PE2(2, pass_a[187], pass_b[194], pass_cs[7][21], pass_a[211], pass_b[218], pass_cs[7][22], final_R[226]);
	PE2_tail(pass_a[164], pass_b[170], pass_cs[6][23], pass_b[195], final_R[202]);

	PE1(13, out_a[18], pass_b[405], out_a[19], PE1_pass_cs[19], out_cs_Q[19], final_R[437]);
	PE2(12, pass_a[375], pass_b[393], pass_cs[18][0], pass_a[388], pass_b[406], pass_cs[18][1], final_R[425]);
	PE2(11, pass_a[363], pass_b[380], pass_cs[17][2], pass_a[377], pass_b[394], pass_cs[17][3], final_R[412]);
	PE2(10, pass_a[350], pass_b[366], pass_cs[16][4], pass_a[365], pass_b[381], pass_cs[16][5], final_R[398]);
	PE2(9, pass_a[336], pass_b[351], pass_cs[15][6], pass_a[352], pass_b[367], pass_cs[15][7], final_R[383]);
	PE2(8, pass_a[321], pass_b[335], pass_cs[14][8], pass_a[338], pass_b[352], pass_cs[14][9], final_R[367]);
	PE2(7, pass_a[305], pass_b[318], pass_cs[13][10], pass_a[323], pass_b[336], pass_cs[13][11], final_R[350]);
	PE2(6, pass_a[288], pass_b[300], pass_cs[12][12], pass_a[307], pass_b[319], pass_cs[12][13], final_R[332]);
	PE2(5, pass_a[270], pass_b[281], pass_cs[11][14], pass_a[290], pass_b[301], pass_cs[11][15], final_R[313]);
	PE2(4, pass_a[251], pass_b[261], pass_cs[10][16], pass_a[272], pass_b[282], pass_cs[10][17], final_R[293]);
	PE2(3, pass_a[231], pass_b[240], pass_cs[9][18], pass_a[253], pass_b[262], pass_cs[9][19], final_R[272]);
	PE2(2, pass_a[210], pass_b[218], pass_cs[8][20], pass_a[233], pass_b[241], pass_cs[8][21], final_R[250]);
	PE2_tail(pass_a[188], pass_b[195], pass_cs[7][22], pass_b[219], final_R[227]);

	PE3(19,out_cs_Q[19],Q_left[18], Q_right[18],Q_left[19], Q_right[19], final_Q_left[19], final_Q_right[19]);
	PE2(12, pass_a[387], pass_b[406], PE1_pass_cs[19], pass_a[399], pass_b[418], pass_cs[19][0], final_R[438]);
	PE2(11, pass_a[376], pass_b[394], pass_cs[18][1], pass_a[389], pass_b[407], pass_cs[18][2], final_R[426]);
	PE2(10, pass_a[364], pass_b[381], pass_cs[17][3], pass_a[378], pass_b[395], pass_cs[17][4], final_R[413]);
	PE2(9, pass_a[351], pass_b[367], pass_cs[16][5], pass_a[366], pass_b[382], pass_cs[16][6], final_R[399]);
	PE2(8, pass_a[337], pass_b[352], pass_cs[15][7], pass_a[353], pass_b[368], pass_cs[15][8], final_R[384]);
	PE2(7, pass_a[322], pass_b[336], pass_cs[14][9], pass_a[339], pass_b[353], pass_cs[14][10], final_R[368]);
	PE2(6, pass_a[306], pass_b[319], pass_cs[13][11], pass_a[324], pass_b[337], pass_cs[13][12], final_R[351]);
	PE2(5, pass_a[289], pass_b[301], pass_cs[12][13], pass_a[308], pass_b[320], pass_cs[12][14], final_R[333]);
	PE2(4, pass_a[271], pass_b[282], pass_cs[11][15], pass_a[291], pass_b[302], pass_cs[11][16], final_R[314]);
	PE2(3, pass_a[252], pass_b[262], pass_cs[10][17], pass_a[273], pass_b[283], pass_cs[10][18], final_R[294]);
	PE2(2, pass_a[232], pass_b[241], pass_cs[9][19], pass_a[254], pass_b[263], pass_cs[9][20], final_R[273]);
	PE2_tail(pass_a[211], pass_b[219], pass_cs[8][21], pass_b[242], final_R[251]);

	PE1(12, out_a[19], pass_b[418], out_a[20], PE1_pass_cs[20], out_cs_Q[20], final_R[450]);
	PE2(11, pass_a[388], pass_b[407], pass_cs[19][0], pass_a[400], pass_b[419], pass_cs[19][1], final_R[439]);
	PE2(10, pass_a[377], pass_b[395], pass_cs[18][2], pass_a[390], pass_b[408], pass_cs[18][3], final_R[427]);
	PE2(9, pass_a[365], pass_b[382], pass_cs[17][4], pass_a[379], pass_b[396], pass_cs[17][5], final_R[414]);
	PE2(8, pass_a[352], pass_b[368], pass_cs[16][6], pass_a[367], pass_b[383], pass_cs[16][7], final_R[400]);
	PE2(7, pass_a[338], pass_b[353], pass_cs[15][8], pass_a[354], pass_b[369], pass_cs[15][9], final_R[385]);
	PE2(6, pass_a[323], pass_b[337], pass_cs[14][10], pass_a[340], pass_b[354], pass_cs[14][11], final_R[369]);
	PE2(5, pass_a[307], pass_b[320], pass_cs[13][12], pass_a[325], pass_b[338], pass_cs[13][13], final_R[352]);
	PE2(4, pass_a[290], pass_b[302], pass_cs[12][14], pass_a[309], pass_b[321], pass_cs[12][15], final_R[334]);
	PE2(3, pass_a[272], pass_b[283], pass_cs[11][16], pass_a[292], pass_b[303], pass_cs[11][17], final_R[315]);
	PE2(2, pass_a[253], pass_b[263], pass_cs[10][18], pass_a[274], pass_b[284], pass_cs[10][19], final_R[295]);
	PE2_tail(pass_a[233], pass_b[242], pass_cs[9][20], pass_b[264], final_R[274]);

	PE3(20,out_cs_Q[20],Q_left[19], Q_right[19],Q_left[20], Q_right[20], final_Q_left[20], final_Q_right[20]);
	PE2(11, pass_a[399], pass_b[419], PE1_pass_cs[20], pass_a[410], pass_b[430], pass_cs[20][0], final_R[451]);
	PE2(10, pass_a[389], pass_b[408], pass_cs[19][1], pass_a[401], pass_b[420], pass_cs[19][2], final_R[440]);
	PE2(9, pass_a[378], pass_b[396], pass_cs[18][3], pass_a[391], pass_b[409], pass_cs[18][4], final_R[428]);
	PE2(8, pass_a[366], pass_b[383], pass_cs[17][5], pass_a[380], pass_b[397], pass_cs[17][6], final_R[415]);
	PE2(7, pass_a[353], pass_b[369], pass_cs[16][7], pass_a[368], pass_b[384], pass_cs[16][8], final_R[401]);
	PE2(6, pass_a[339], pass_b[354], pass_cs[15][9], pass_a[355], pass_b[370], pass_cs[15][10], final_R[386]);
	PE2(5, pass_a[324], pass_b[338], pass_cs[14][11], pass_a[341], pass_b[355], pass_cs[14][12], final_R[370]);
	PE2(4, pass_a[308], pass_b[321], pass_cs[13][13], pass_a[326], pass_b[339], pass_cs[13][14], final_R[353]);
	PE2(3, pass_a[291], pass_b[303], pass_cs[12][15], pass_a[310], pass_b[322], pass_cs[12][16], final_R[335]);
	PE2(2, pass_a[273], pass_b[284], pass_cs[11][17], pass_a[293], pass_b[304], pass_cs[11][18], final_R[316]);
	PE2_tail(pass_a[254], pass_b[264], pass_cs[10][19], pass_b[285], final_R[296]);

	PE1(11, out_a[20], pass_b[430], out_a[21], PE1_pass_cs[21], out_cs_Q[21], final_R[462]);
	PE2(10, pass_a[400], pass_b[420], pass_cs[20][0], pass_a[411], pass_b[431], pass_cs[20][1], final_R[452]);
	PE2(9, pass_a[390], pass_b[409], pass_cs[19][2], pass_a[402], pass_b[421], pass_cs[19][3], final_R[441]);
	PE2(8, pass_a[379], pass_b[397], pass_cs[18][4], pass_a[392], pass_b[410], pass_cs[18][5], final_R[429]);
	PE2(7, pass_a[367], pass_b[384], pass_cs[17][6], pass_a[381], pass_b[398], pass_cs[17][7], final_R[416]);
	PE2(6, pass_a[354], pass_b[370], pass_cs[16][8], pass_a[369], pass_b[385], pass_cs[16][9], final_R[402]);
	PE2(5, pass_a[340], pass_b[355], pass_cs[15][10], pass_a[356], pass_b[371], pass_cs[15][11], final_R[387]);
	PE2(4, pass_a[325], pass_b[339], pass_cs[14][12], pass_a[342], pass_b[356], pass_cs[14][13], final_R[371]);
	PE2(3, pass_a[309], pass_b[322], pass_cs[13][14], pass_a[327], pass_b[340], pass_cs[13][15], final_R[354]);
	PE2(2, pass_a[292], pass_b[304], pass_cs[12][16], pass_a[311], pass_b[323], pass_cs[12][17], final_R[336]);
	PE2_tail(pass_a[274], pass_b[285], pass_cs[11][18], pass_b[305], final_R[317]);

	PE3(21,out_cs_Q[21],Q_left[20], Q_right[20],Q_left[21], Q_right[21], final_Q_left[21], final_Q_right[21]);
	PE2(10, pass_a[410], pass_b[431], PE1_pass_cs[21], pass_a[420], pass_b[441], pass_cs[21][0], final_R[463]);
	PE2(9, pass_a[401], pass_b[421], pass_cs[20][1], pass_a[412], pass_b[432], pass_cs[20][2], final_R[453]);
	PE2(8, pass_a[391], pass_b[410], pass_cs[19][3], pass_a[403], pass_b[422], pass_cs[19][4], final_R[442]);
	PE2(7, pass_a[380], pass_b[398], pass_cs[18][5], pass_a[393], pass_b[411], pass_cs[18][6], final_R[430]);
	PE2(6, pass_a[368], pass_b[385], pass_cs[17][7], pass_a[382], pass_b[399], pass_cs[17][8], final_R[417]);
	PE2(5, pass_a[355], pass_b[371], pass_cs[16][9], pass_a[370], pass_b[386], pass_cs[16][10], final_R[403]);
	PE2(4, pass_a[341], pass_b[356], pass_cs[15][11], pass_a[357], pass_b[372], pass_cs[15][12], final_R[388]);
	PE2(3, pass_a[326], pass_b[340], pass_cs[14][13], pass_a[343], pass_b[357], pass_cs[14][14], final_R[372]);
	PE2(2, pass_a[310], pass_b[323], pass_cs[13][15], pass_a[328], pass_b[341], pass_cs[13][16], final_R[355]);
	PE2_tail(pass_a[293], pass_b[305], pass_cs[12][17], pass_b[324], final_R[337]);

	PE1(10, out_a[21], pass_b[441], out_a[22], PE1_pass_cs[22], out_cs_Q[22], final_R[473]);
	PE2(9, pass_a[411], pass_b[432], pass_cs[21][0], pass_a[421], pass_b[442], pass_cs[21][1], final_R[464]);
	PE2(8, pass_a[402], pass_b[422], pass_cs[20][2], pass_a[413], pass_b[433], pass_cs[20][3], final_R[454]);
	PE2(7, pass_a[392], pass_b[411], pass_cs[19][4], pass_a[404], pass_b[423], pass_cs[19][5], final_R[443]);
	PE2(6, pass_a[381], pass_b[399], pass_cs[18][6], pass_a[394], pass_b[412], pass_cs[18][7], final_R[431]);
	PE2(5, pass_a[369], pass_b[386], pass_cs[17][8], pass_a[383], pass_b[400], pass_cs[17][9], final_R[418]);
	PE2(4, pass_a[356], pass_b[372], pass_cs[16][10], pass_a[371], pass_b[387], pass_cs[16][11], final_R[404]);
	PE2(3, pass_a[342], pass_b[357], pass_cs[15][12], pass_a[358], pass_b[373], pass_cs[15][13], final_R[389]);
	PE2(2, pass_a[327], pass_b[341], pass_cs[14][14], pass_a[344], pass_b[358], pass_cs[14][15], final_R[373]);
	PE2_tail(pass_a[311], pass_b[324], pass_cs[13][16], pass_b[342], final_R[356]);

	PE3(22,out_cs_Q[22],Q_left[21], Q_right[21],Q_left[22], Q_right[22], final_Q_left[22], final_Q_right[22]);
	PE2(9, pass_a[420], pass_b[442], PE1_pass_cs[22], pass_a[429], pass_b[451], pass_cs[22][0], final_R[474]);
	PE2(8, pass_a[412], pass_b[433], pass_cs[21][1], pass_a[422], pass_b[443], pass_cs[21][2], final_R[465]);
	PE2(7, pass_a[403], pass_b[423], pass_cs[20][3], pass_a[414], pass_b[434], pass_cs[20][4], final_R[455]);
	PE2(6, pass_a[393], pass_b[412], pass_cs[19][5], pass_a[405], pass_b[424], pass_cs[19][6], final_R[444]);
	PE2(5, pass_a[382], pass_b[400], pass_cs[18][7], pass_a[395], pass_b[413], pass_cs[18][8], final_R[432]);
	PE2(4, pass_a[370], pass_b[387], pass_cs[17][9], pass_a[384], pass_b[401], pass_cs[17][10], final_R[419]);
	PE2(3, pass_a[357], pass_b[373], pass_cs[16][11], pass_a[372], pass_b[388], pass_cs[16][12], final_R[405]);
	PE2(2, pass_a[343], pass_b[358], pass_cs[15][13], pass_a[359], pass_b[374], pass_cs[15][14], final_R[390]);
	PE2_tail(pass_a[328], pass_b[342], pass_cs[14][15], pass_b[359], final_R[374]);

	PE1(9, out_a[22], pass_b[451], out_a[23], PE1_pass_cs[23], out_cs_Q[23], final_R[483]);
	PE2(8, pass_a[421], pass_b[443], pass_cs[22][0], pass_a[430], pass_b[452], pass_cs[22][1], final_R[475]);
	PE2(7, pass_a[413], pass_b[434], pass_cs[21][2], pass_a[423], pass_b[444], pass_cs[21][3], final_R[466]);
	PE2(6, pass_a[404], pass_b[424], pass_cs[20][4], pass_a[415], pass_b[435], pass_cs[20][5], final_R[456]);
	PE2(5, pass_a[394], pass_b[413], pass_cs[19][6], pass_a[406], pass_b[425], pass_cs[19][7], final_R[445]);
	PE2(4, pass_a[383], pass_b[401], pass_cs[18][8], pass_a[396], pass_b[414], pass_cs[18][9], final_R[433]);
	PE2(3, pass_a[371], pass_b[388], pass_cs[17][10], pass_a[385], pass_b[402], pass_cs[17][11], final_R[420]);
	PE2(2, pass_a[358], pass_b[374], pass_cs[16][12], pass_a[373], pass_b[389], pass_cs[16][13], final_R[406]);
	PE2_tail(pass_a[344], pass_b[359], pass_cs[15][14], pass_b[375], final_R[391]);

	PE3(23,out_cs_Q[23],Q_left[22], Q_right[22],Q_left[23], Q_right[23], final_Q_left[23], final_Q_right[23]);
	PE2(8, pass_a[429], pass_b[452], PE1_pass_cs[23], pass_a[437], pass_b[460], pass_cs[23][0], final_R[484]);
	PE2(7, pass_a[422], pass_b[444], pass_cs[22][1], pass_a[431], pass_b[453], pass_cs[22][2], final_R[476]);
	PE2(6, pass_a[414], pass_b[435], pass_cs[21][3], pass_a[424], pass_b[445], pass_cs[21][4], final_R[467]);
	PE2(5, pass_a[405], pass_b[425], pass_cs[20][5], pass_a[416], pass_b[436], pass_cs[20][6], final_R[457]);
	PE2(4, pass_a[395], pass_b[414], pass_cs[19][7], pass_a[407], pass_b[426], pass_cs[19][8], final_R[446]);
	PE2(3, pass_a[384], pass_b[402], pass_cs[18][9], pass_a[397], pass_b[415], pass_cs[18][10], final_R[434]);
	PE2(2, pass_a[372], pass_b[389], pass_cs[17][11], pass_a[386], pass_b[403], pass_cs[17][12], final_R[421]);
	PE2_tail(pass_a[359], pass_b[375], pass_cs[16][13], pass_b[390], final_R[407]);

	PE1(8, out_a[23], pass_b[460], out_a[24], PE1_pass_cs[24], out_cs_Q[24], final_R[492]);
	PE2(7, pass_a[430], pass_b[453], pass_cs[23][0], pass_a[438], pass_b[461], pass_cs[23][1], final_R[485]);
	PE2(6, pass_a[423], pass_b[445], pass_cs[22][2], pass_a[432], pass_b[454], pass_cs[22][3], final_R[477]);
	PE2(5, pass_a[415], pass_b[436], pass_cs[21][4], pass_a[425], pass_b[446], pass_cs[21][5], final_R[468]);
	PE2(4, pass_a[406], pass_b[426], pass_cs[20][6], pass_a[417], pass_b[437], pass_cs[20][7], final_R[458]);
	PE2(3, pass_a[396], pass_b[415], pass_cs[19][8], pass_a[408], pass_b[427], pass_cs[19][9], final_R[447]);
	PE2(2, pass_a[385], pass_b[403], pass_cs[18][10], pass_a[398], pass_b[416], pass_cs[18][11], final_R[435]);
	PE2_tail(pass_a[373], pass_b[390], pass_cs[17][12], pass_b[404], final_R[422]);

	PE3(24,out_cs_Q[24],Q_left[23], Q_right[23],Q_left[24], Q_right[24], final_Q_left[24], final_Q_right[24]);
	PE2(7, pass_a[437], pass_b[461], PE1_pass_cs[24], pass_a[444], pass_b[468], pass_cs[24][0], final_R[493]);
	PE2(6, pass_a[431], pass_b[454], pass_cs[23][1], pass_a[439], pass_b[462], pass_cs[23][2], final_R[486]);
	PE2(5, pass_a[424], pass_b[446], pass_cs[22][3], pass_a[433], pass_b[455], pass_cs[22][4], final_R[478]);
	PE2(4, pass_a[416], pass_b[437], pass_cs[21][5], pass_a[426], pass_b[447], pass_cs[21][6], final_R[469]);
	PE2(3, pass_a[407], pass_b[427], pass_cs[20][7], pass_a[418], pass_b[438], pass_cs[20][8], final_R[459]);
	PE2(2, pass_a[397], pass_b[416], pass_cs[19][9], pass_a[409], pass_b[428], pass_cs[19][10], final_R[448]);
	PE2_tail(pass_a[386], pass_b[404], pass_cs[18][11], pass_b[417], final_R[436]);

	PE1(7, out_a[24], pass_b[468], out_a[25], PE1_pass_cs[25], out_cs_Q[25], final_R[500]);
	PE2(6, pass_a[438], pass_b[462], pass_cs[24][0], pass_a[445], pass_b[469], pass_cs[24][1], final_R[494]);
	PE2(5, pass_a[432], pass_b[455], pass_cs[23][2], pass_a[440], pass_b[463], pass_cs[23][3], final_R[487]);
	PE2(4, pass_a[425], pass_b[447], pass_cs[22][4], pass_a[434], pass_b[456], pass_cs[22][5], final_R[479]);
	PE2(3, pass_a[417], pass_b[438], pass_cs[21][6], pass_a[427], pass_b[448], pass_cs[21][7], final_R[470]);
	PE2(2, pass_a[408], pass_b[428], pass_cs[20][8], pass_a[419], pass_b[439], pass_cs[20][9], final_R[460]);
	PE2_tail(pass_a[398], pass_b[417], pass_cs[19][10], pass_b[429], final_R[449]);

	PE3(25,out_cs_Q[25],Q_left[24], Q_right[24],Q_left[25], Q_right[25], final_Q_left[25], final_Q_right[25]);
	PE2(6, pass_a[444], pass_b[469], PE1_pass_cs[25], pass_a[450], pass_b[475], pass_cs[25][0], final_R[501]);
	PE2(5, pass_a[439], pass_b[463], pass_cs[24][1], pass_a[446], pass_b[470], pass_cs[24][2], final_R[495]);
	PE2(4, pass_a[433], pass_b[456], pass_cs[23][3], pass_a[441], pass_b[464], pass_cs[23][4], final_R[488]);
	PE2(3, pass_a[426], pass_b[448], pass_cs[22][5], pass_a[435], pass_b[457], pass_cs[22][6], final_R[480]);
	PE2(2, pass_a[418], pass_b[439], pass_cs[21][7], pass_a[428], pass_b[449], pass_cs[21][8], final_R[471]);
	PE2_tail(pass_a[409], pass_b[429], pass_cs[20][9], pass_b[440], final_R[461]);

	PE1(6, out_a[25], pass_b[475], out_a[26], PE1_pass_cs[26], out_cs_Q[26], final_R[507]);
	PE2(5, pass_a[445], pass_b[470], pass_cs[25][0], pass_a[451], pass_b[476], pass_cs[25][1], final_R[502]);
	PE2(4, pass_a[440], pass_b[464], pass_cs[24][2], pass_a[447], pass_b[471], pass_cs[24][3], final_R[496]);
	PE2(3, pass_a[434], pass_b[457], pass_cs[23][4], pass_a[442], pass_b[465], pass_cs[23][5], final_R[489]);
	PE2(2, pass_a[427], pass_b[449], pass_cs[22][6], pass_a[436], pass_b[458], pass_cs[22][7], final_R[481]);
	PE2_tail(pass_a[419], pass_b[440], pass_cs[21][8], pass_b[450], final_R[472]);

	PE3(26,out_cs_Q[26],Q_left[25], Q_right[25],Q_left[26], Q_right[26], final_Q_left[26], final_Q_right[26]);
	PE2(5, pass_a[450], pass_b[476], PE1_pass_cs[26], pass_a[455], pass_b[481], pass_cs[26][0], final_R[508]);
	PE2(4, pass_a[446], pass_b[471], pass_cs[25][1], pass_a[452], pass_b[477], pass_cs[25][2], final_R[503]);
	PE2(3, pass_a[441], pass_b[465], pass_cs[24][3], pass_a[448], pass_b[472], pass_cs[24][4], final_R[497]);
	PE2(2, pass_a[435], pass_b[458], pass_cs[23][5], pass_a[443], pass_b[466], pass_cs[23][6], final_R[490]);
	PE2_tail(pass_a[428], pass_b[450], pass_cs[22][7], pass_b[459], final_R[482]);

	PE1(5, out_a[26], pass_b[481], out_a[27], PE1_pass_cs[27], out_cs_Q[27], final_R[513]);
	PE2(4, pass_a[451], pass_b[477], pass_cs[26][0], pass_a[456], pass_b[482], pass_cs[26][1], final_R[509]);
	PE2(3, pass_a[447], pass_b[472], pass_cs[25][2], pass_a[453], pass_b[478], pass_cs[25][3], final_R[504]);
	PE2(2, pass_a[442], pass_b[466], pass_cs[24][4], pass_a[449], pass_b[473], pass_cs[24][5], final_R[498]);
	PE2_tail(pass_a[436], pass_b[459], pass_cs[23][6], pass_b[467], final_R[491]);

	PE3(27,out_cs_Q[27],Q_left[26], Q_right[26],Q_left[27], Q_right[27], final_Q_left[27], final_Q_right[27]);
	PE2(4, pass_a[455], pass_b[482], PE1_pass_cs[27], pass_a[459], pass_b[486], pass_cs[27][0], final_R[514]);
	PE2(3, pass_a[452], pass_b[478], pass_cs[26][1], pass_a[457], pass_b[483], pass_cs[26][2], final_R[510]);
	PE2(2, pass_a[448], pass_b[473], pass_cs[25][3], pass_a[454], pass_b[479], pass_cs[25][4], final_R[505]);
	PE2_tail(pass_a[443], pass_b[467], pass_cs[24][5], pass_b[474], final_R[499]);

	PE1(4, out_a[27], pass_b[486], out_a[28], PE1_pass_cs[28], out_cs_Q[28], final_R[518]);
	PE2(3, pass_a[456], pass_b[483], pass_cs[27][0], pass_a[460], pass_b[487], pass_cs[27][1], final_R[515]);
	PE2(2, pass_a[453], pass_b[479], pass_cs[26][2], pass_a[458], pass_b[484], pass_cs[26][3], final_R[511]);
	PE2_tail(pass_a[449], pass_b[474], pass_cs[25][4], pass_b[480], final_R[506]);

	PE3(28,out_cs_Q[28],Q_left[27], Q_right[27],Q_left[28], Q_right[28], final_Q_left[28], final_Q_right[28]);
	PE2(3, pass_a[459], pass_b[487], PE1_pass_cs[28], pass_a[462], pass_b[490], pass_cs[28][0], final_R[519]);
	PE2(2, pass_a[457], pass_b[484], pass_cs[27][1], pass_a[461], pass_b[488], pass_cs[27][2], final_R[516]);
	PE2_tail(pass_a[454], pass_b[480], pass_cs[26][3], pass_b[485], final_R[512]);

	PE1(3, out_a[28], pass_b[490], out_a[29], PE1_pass_cs[29], out_cs_Q[29], final_R[522]);
	PE2(2, pass_a[460], pass_b[488], pass_cs[28][0], pass_a[463], pass_b[491], pass_cs[28][1], final_R[520]);
	PE2_tail(pass_a[458], pass_b[485], pass_cs[27][2], pass_b[489], final_R[517]);

	PE3(29,out_cs_Q[29],Q_left[28], Q_right[28],Q_left[29], Q_right[29], final_Q_left[29], final_Q_right[29]);
	PE2(2, pass_a[462], pass_b[491], PE1_pass_cs[29], pass_a[464], pass_b[493], pass_cs[29][0], final_R[523]);
	PE2_tail(pass_a[461], pass_b[489], pass_cs[28][1], pass_b[492], final_R[521]);

	PE1(2, out_a[29], pass_b[493], out_a[30], PE1_pass_cs[30], out_cs_Q[30], final_R[525]);
	PE2_tail(pass_a[463], pass_b[492], pass_cs[29][0], pass_b[494], final_R[524]);

	PE3(30,out_cs_Q[30],Q_left[29], Q_right[29],Q_left[30], Q_right[30], final_Q_left[30], final_Q_right[30]);
	PE2_tail(pass_a[464], pass_b[494], PE1_pass_cs[30], pass_b[495], final_R[526]);

	PE1_tail(out_a[30], pass_b[495], out_cs_Q[31], final_R[527]);

	PE3_tail(31,out_cs_Q[31], Q_left[30], Q_right[30], final_Q_left[31], final_Q_right[31]);

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
  
