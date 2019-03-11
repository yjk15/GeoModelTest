#pragma once
#include <QString>
#include "matrix.h"
#include <vector>
using namespace std;

class MODEL {
public:
	int testType, model;
	
	int endAndReversalType; //���ü���������߷�תʱ��״̬��0-->p; 1-->q; 2-->��Ӧ��
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