#pragma once
#include <QString>
#include "matrix.h"
#include <vector>
#include <cmath>
using namespace std;

class MODEL {
public:
	int testType, model;
	
	double ee; //ģ�⿪ʼʱ�Ŀ�϶��
	int endAndReversalType; //���ü���������߷�תʱ��״̬��0-->p; 1-->q; 2-->��Ӧ��
	double endAndReversalPoint;
	int loop, loopCounter;
	bool direction;

	QString *figureTitle;
	int axisX, axisY;
	double stepLength;
	double *internalParameter;
	double pAtmos;
	vector<MATRIX> *stressPath, *strainPath;

	MATRIX stress, strain;
	MATRIX stressIncrement, strainIncrement;
	vector<vector<double>> saveParameter; //0���϶��e��1-9��alpha��10-18��z, 19-27��alphaInit

public:
	MODEL();
	~MODEL();

	void Simulate();

private:
	void GetStrainIncrementForSpecifiedTestType();
	void Integrator(bool);
	bool isEndingPoint();
	bool isReversalPoint();

	void IntegratorE(bool);

	void IntegratorDMExplicit(bool);
	struct RK4Class {
		MATRIX ds, dAlpha, dz;
		double dee;
	};
	RK4Class RK4(MATRIX stress, MATRIX alpha, double ee, MATRIX z, MATRIX srtain, MATRIX alphaInit);
	void IntegratorDMImplicit(bool);
	double relu(double x);
	double getG(double p, double e);
	double getK(double G);
	double getEc(double pc);
	double getF(MATRIX s, MATRIX alpha, double p);
	MATRIX getN(MATRIX r, MATRIX alpha);
	double getCos3Theta(MATRIX n);
	double getg(double cos3Theta);
	MATRIX getAlphaThetaD(MATRIX n, double p, double e, double g);
	double getB(double cos3Theta, double g);
	double getC(double g);
	double getAd(MATRIX z, MATRIX n);
	double getD(MATRIX n, double Ad, MATRIX alpha, MATRIX alphaThetaD);
	MATRIX getAlphaThetaB(MATRIX n, double p, double e, double g);
	double getH(double e, double p, MATRIX alpha, MATRIX alphaInit, MATRIX n);
	double getKp(MATRIX alphaThetaB, MATRIX alpha, double p, MATRIX n, double h);
	double getL(MATRIX n, double G, MATRIX r, MATRIX de, double depsv, double Kp, double B, double C, double K, double D);
	MATRIX getRAp(double B, double C, MATRIX n);
	MATRIX getdSigma(double G, MATRIX de, double K, double depsv, double L, MATRIX RAp, double D);
	MATRIX getdAlpha(double L, double h, MATRIX alphaThetaB, MATRIX alpha);
	MATRIX getdz(double depsvp, MATRIX n, MATRIX z);
};