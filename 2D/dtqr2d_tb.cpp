#include "input.h"
#include "dtqr2d.h"
using namespace std;

void PrintMatrix(double A[ROWS][COLS])
{
	for (int i = 0; i < ROWS; ++i)
	{
		for (int j = 0; j < COLS; ++j)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void PrintMatrixSquare(double A[ROWS][ROWS])
{
	for (int i = 0; i < ROWS; ++i)
	{
		for (int j = 0; j < ROWS; ++j)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void MM(double A[ROWS][ROWS], double B[ROWS][COLS], double C[ROWS][COLS])
{
	for (int i = 0; i < ROWS; ++i)
	{
		for (int j = 0; j < COLS; ++j)
		{
			double res = 0;
			for (int k = 0; k < COLS; ++k)
			{
				res += double(A[i][k]) * double(B[k][j]);
			}
			C[i][j] = res;
		}
	}

	cout << endl;
}

void MakeEye(double A[ROWS][ROWS])
{
	for (int i = 0; i < ROWS; ++i)
	{
		for (int j = 0; j < ROWS; ++j)
		{
			if (i == j)
				A[i][j] = 1;
			else
				A[i][j] = 0;
		}
	}
}

void stream2matrix(MATRIX_T A1[LEN], MATRIX_T A2[LEN], double Target_A[ROWS][COLS])
{
	int count = 0;
	for (int i = 0; i < COLS; ++i)
	{
		for (int j = i; j < COLS; ++j)
		{
			Target_A[i][j] = A1[count];
			count++;
		}
	}

	count = 0;
	for (int i = COLS; i < ROWS; ++i)
	{
		for (int j = i - COLS; j < COLS; ++j)
		{
			Target_A[i][j] = A2[count];
			count++;
		}
	}
}

void TestEachBench(MATRIX_T A1[LEN], MATRIX_T A2[LEN])
{
	int qr_GR_success = 0;

	double Target_A[ROWS][COLS] = {0};
	double Q[ROWS][ROWS] = {0};

	MakeEye(Q);
	fixed_cs Q_L[Q_LEN];
	fixed_cs Q_R[Q_LEN];
	MATRIX_T R[LEN] = {0};

	double R_re[ROWS][COLS] = {0};
	double A_re[ROWS][COLS] = {0};

	stream<MATRIX_T> A1_s;
	stream<MATRIX_T> A2_s;
	stream<fixed_cs> Q_L_s;
	stream<fixed_cs> Q_R_s;
	stream<MATRIX_T> R_s;


	cout << "Target Matrix" << endl;
	stream2matrix(A1, A2, Target_A);
	PrintMatrix(Target_A);

	for (int i=0;i<LEN;i++){
		A1_s.write(A1[i]);
		A2_s.write(A2[i]);
	}
	qr_GR_success = top(A1_s, A2_s, Q_L_s, Q_R_s, R_s);
	for (int i=0;i<Q_LEN;i++){
		Q_L[i]=Q_L_s.read();
		Q_R[i]=Q_R_s.read();
	}
	int count = 0;
	for (int i = 0; i < COLS; ++i)
	{
		for (int j = i; j < COLS; ++j)
		{
			R_re[i][j] = R_s.read();
			count += 1;
		}
	}
	int counts = 0;
	for (int i = 0; i < COLS; ++i)
	{
		for (int j = 0; j < i + 1; ++j)
		{
			Q[j][i] = Q_L[counts];
			Q[COLS - i - 1 + j][COLS * 2 - 1 - i] = Q_R[counts];
			counts++;
		}
		for (int j = 0; j < i + 1; ++j)
		{
			Q[j + COLS][i] = Q_L[counts];
			Q[COLS * 2 - i - 1 + j][COLS * 2 - 1 - i] = Q_R[counts];
			counts++;
		}
	}
	cout << "Matrix R: " << endl;
	PrintMatrix(R_re);

	cout << "Matrix Q: " << endl;
	PrintMatrixSquare(Q);

	MM(Q, R_re, A_re);
	cout << "Results (Q * R): " << endl;
	PrintMatrix(A_re);

	double sum_err_rate = 0.0;
	double max_err = 0.0;
	for(int i = 0; i < ROWS; ++i)
	{
		for(int j = 0; j < COLS; ++j)
		{	
			double err = hls::abs(A_re[i][j] - Target_A[i][j]);
			double err_rate = Target_A[i][j] == 0 ? 0 : err / Target_A[i][j];
			sum_err_rate += err_rate;
			if (err > max_err)
				max_err = err;
		}
	}
	cout<<sum_err_rate<<endl<<max_err<<endl;


}
int main()
{

	cout << endl;
	cout << "First Test Bench" << endl;
	MATRIX_T A1[LEN];
	MATRIX_T A2[LEN];
	int* A1_in;
	int* A2_in;

	switch (COLS){
		case 4:
			A1_in = A1_4;
			A2_in = A2_4;
			break;
		case 8:
			A1_in = A1_8;
			A2_in = A2_8;
			break;
		case 12:
			A1_in = A1_12;
			A2_in = A2_12;
			break;
		case 16:
			A1_in = A1_16;
			A2_in = A2_16;
			break;
		case 24:
			A1_in = A1_24;
			A2_in = A2_24;
			break;
		default: 
        	cout<<"error size"<<endl;
        	break;
	}

	for (int i = 0; i < LEN; ++i){
		A1[i] = A1_in[i];
		A2[i] = A2_in[i];
	}
	
	TestEachBench(A1, A2);
	return 0;
}
