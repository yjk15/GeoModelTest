#pragma once
#include <QString>
#include "matrix.h"
#include <vector>
using namespace std;

class MODEL {
public:
	int testType, model;
	
	int endAndReversalType; //设置计算结束或者反转时的状态，0-->p; 1-->q; 2-->体应变
	double endAndReversalPoint;
	int loop, loopCounter;
	bool direction;

	QString *figureTitle;
	int axisX, axisY;
	double stepLength;
	double *internalParameter;
	vector<MATRIX> *stressPath, *strainPath;

	MATRIX stress, strain;
	MATRIX stressIncrement, strainIncrement;

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
};