#include "matrix.h"

MATRIX::MATRIX() {
	this->matrix = new double[9];
	for (int i = 0; i < 9; i++)
		this->matrix[i] = 0;
}

MATRIX::MATRIX(double a1, double a2, double a3) {
	this->matrix = new double[9];
	for (int i = 0; i < 9; i++)
		this->matrix[i] = 0;
	this->matrix[0] = a1;
	this->matrix[4] = a2;
	this->matrix[8] = a3;
}

MATRIX::MATRIX(double a[9]) {
	this->matrix = new double[9];
	for (int i = 0; i < 9; i++)
		this->matrix[i] = a[i];
}

MATRIX::MATRIX(const MATRIX& M) {
	this->matrix = new double[9];
	for (int i = 0; i < 9; i++)
		this->matrix[i] = M.matrix[i];
}

MATRIX::~MATRIX() {
	delete[] this->matrix;
}

double& MATRIX::operator()(int i, int j) {
	if (i > 2 || i < 0 || j > 2 || j < 0)
		return this->matrix[0];
	return this->matrix[i * 3 + j];
}

MATRIX& MATRIX::operator=(const MATRIX& M) {
	for (int i = 0; i < 9; i++)
		this->matrix[i] = M.matrix[i];

	return *this;
}

void MATRIX::clear() {
	for (int i = 0; i < 9; i++)
		this->matrix[i] = 0;
}

MATRIX operator+(MATRIX M, MATRIX N) {
	MATRIX tmp = M;
	for (int i = 0; i < 9; i++)
		tmp.matrix[i] += N.matrix[i];

	return tmp;
}

MATRIX operator-(MATRIX M, MATRIX N) {
	MATRIX tmp = M;
	for (int i = 0; i < 9; i++)
		tmp.matrix[i] -= N.matrix[i];

	return tmp;
}

MATRIX operator * (MATRIX M, MATRIX N) {
	MATRIX tmp;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				tmp(i, k) += M(i, j) * N(j, k);

	return tmp;
}

MATRIX operator * (MATRIX M, double a) {
	MATRIX tmp = M;
	for (int i = 0; i < 9; i++)
		tmp.matrix[i] *= a;

	return tmp;
}

MATRIX operator * (double a, MATRIX M) {
	MATRIX tmp = M;
	for (int i = 0; i < 9; i++)
		tmp.matrix[i] *= a;

	return tmp;
}

MATRIX operator / (MATRIX M, double a) {
	MATRIX tmp = M;
	for (int i = 0; i < 9; i++)
		tmp.matrix[i] /= a;

	return tmp;
}

double operator % (MATRIX M, MATRIX N) {
	double a = 0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			a += M(i, j) * N(j, i);

	return a;
}

void print(MATRIX M) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << std::setw(6) << M(i, j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

double tr(MATRIX M) {
	return M(0, 0) + M(1, 1) + M(2, 2);
}