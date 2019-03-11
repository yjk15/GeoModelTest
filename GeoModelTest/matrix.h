#pragma once
#include <iostream>
#include <iomanip>

class MATRIX {
public:
	double *matrix;

public:
	MATRIX();
	MATRIX(double a1, double a2, double a3);
	MATRIX(double a[9]);
	MATRIX(const MATRIX& M);
	~MATRIX();

	double& operator () (int i, int j);
	MATRIX& operator = (const MATRIX& M);
	void clear(); //������Ԫ�ر�Ϊ0
};

MATRIX operator + (MATRIX& M, MATRIX& N);
MATRIX operator - (MATRIX& M, MATRIX& N);
MATRIX operator * (MATRIX& M, MATRIX& N);
MATRIX operator * (MATRIX& M, double a);
MATRIX operator * (double a, MATRIX& M);
MATRIX operator / (MATRIX& M, double a);
double operator % (MATRIX& M, MATRIX& N); //�á�%������ʾ��:����������������Ԫ�ض�λ���֮��
void print(  MATRIX& M);
double tr(MATRIX& M);