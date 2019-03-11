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
	void clear(); //将所有元素变为0
};

MATRIX operator + (MATRIX& M, MATRIX& N);
MATRIX operator - (MATRIX& M, MATRIX& N);
MATRIX operator * (MATRIX& M, MATRIX& N);
MATRIX operator * (MATRIX& M, double a);
MATRIX operator * (double a, MATRIX& M);
MATRIX operator / (MATRIX& M, double a);
double operator % (MATRIX& M, MATRIX& N); //用“%”来表示“:”，即两矩阵所有元素对位相乘之和
void print(  MATRIX& M);
double tr(MATRIX& M);