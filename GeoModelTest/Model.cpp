#include "Model.h"

MODEL::MODEL() {
	testType = 0;
	model = 0;

	endAndReversalType = 0;
	endAndReversalPoint = 0;
	loop = 0;
	loopCounter = 0;
	direction = false;

	figureTitle = new QString();
	axisX = 0;
	axisY = 0;

	stepLength = 1e-5;

	internalParameter = new double[30];
	for (int i = 0; i < 30; i++)
		internalParameter[i] = 0;

	stressPath = new vector<MATRIX>;
	strainPath = new vector<MATRIX>;
}

MODEL::~MODEL() {
	delete figureTitle;
	delete[] internalParameter;
	delete stressPath;
	delete strainPath;
}

void MODEL::Simulate() {
	while (!isEndingPoint()) {
		stressPath->push_back(stress);
		strainPath->push_back(strain);
		GetStrainIncrementForSpecifiedTestType();
		Integrator(true);
		strain = strain + strainIncrement;
		stress = stress + stressIncrement;
	}
	stressPath->push_back(stress);
	strainPath->push_back(strain);
}

void MODEL::GetStrainIncrementForSpecifiedTestType() {
	double stressTolerance = 0.3;

	switch (testType) {
	case 0:
		strainIncrement.clear();
		strainIncrement(0, 0) = stepLength / 2;
		strainIncrement(1, 1) = stepLength / 2;
		strainIncrement(2, 2) = -stepLength;
		break;

	case 1:
		strainIncrement.clear();
		strainIncrement(0, 0) = -stepLength / 2;
		strainIncrement(1, 1) = -stepLength / 2;
		strainIncrement(2, 2) = stepLength;
		break;

	case 2:
		strainIncrement.clear();
		if (direction) {
			strainIncrement(0, 0) = stepLength / 2;
			strainIncrement(1, 1) = stepLength / 2;
			strainIncrement(2, 2) = -stepLength;
		}
		else {
			strainIncrement(0, 0) = -stepLength / 2;
			strainIncrement(1, 1) = -stepLength / 2;
			strainIncrement(2, 2) = stepLength;
		}
		break;

	case 3:
		strainIncrement.clear();
		strainIncrement(0, 0) = 0;
		strainIncrement(1, 1) = 0;
		strainIncrement(2, 2) = -stepLength;
		do {
			Integrator(false);
			if (stressIncrement(0, 0) < 0) 
				strainIncrement(0, 0) += stepLength / 100;
			else
				strainIncrement(0, 0) -= stepLength / 100;
			if (stressIncrement(1, 1) < 0)
				strainIncrement(1, 1) += stepLength / 100;
			else
				strainIncrement(1, 1) -= stepLength / 100;
		} while (abs(stressIncrement(0, 0)) > 0 || abs(stressIncrement(1, 1)) > 0);
		break;

	case 4:
		strainIncrement.clear();
		strainIncrement(0, 0) = 0;
		strainIncrement(1, 1) = 0;
		strainIncrement(2, 2) = stepLength;
		do {
			Integrator(false);
			if (stressIncrement(0, 0) < 0)
				strainIncrement(0, 0) += stepLength / 100;
			else
				strainIncrement(0, 0) -= stepLength / 100;
			if (stressIncrement(1, 1) < 0)
				strainIncrement(1, 1) += stepLength / 100;
			else
				strainIncrement(1, 1) -= stepLength / 100;
		} while (abs(stressIncrement(0, 0)) > stressTolerance || abs(stressIncrement(1, 1)) > stressTolerance);
		break;

	case 5:
		strainIncrement.clear();
		strainIncrement(0, 0) = 0;
		strainIncrement(1, 1) = 0;
		if (direction)
			strainIncrement(2, 2) = -stepLength;
		else
			strainIncrement(2, 2) = stepLength;
		do {
			Integrator(false);
			if (stressIncrement(0, 0) < 0)
				strainIncrement(0, 0) += stepLength / 100;
			else
				strainIncrement(0, 0) -= stepLength / 100;
			if (stressIncrement(1, 1) < 0)
				strainIncrement(1, 1) += stepLength / 100;
			else
				strainIncrement(1, 1) -= stepLength / 100;
		} while (abs(stressIncrement(0, 0)) > stressTolerance || abs(stressIncrement(1, 1)) > stressTolerance);
		break;

	default:
		strainIncrement.clear();
		break;
	}
}

void MODEL::Integrator(bool updateFlag) {
	switch (model) {
	case 0:
		IntegratorE(updateFlag);
		break;
	default:
		stressIncrement.clear();
		break;
	}
}

void MODEL::IntegratorE(bool updateFlag) {
	double E = internalParameter[0];
	double v = internalParameter[1];
	double mu2 = E / (1.0 + v);
	double lam = v * mu2 / (1.0 - 2.0*v);
	double mu = 0.50*mu2;

	mu2 += lam;

	stressIncrement(0, 0) = mu2 * strainIncrement(0, 0) + lam * (strainIncrement(1, 1) + strainIncrement(2, 2));
	stressIncrement(1, 1) = mu2 * strainIncrement(1, 1) + lam * (strainIncrement(0, 0) + strainIncrement(2, 2));
	stressIncrement(2, 2) = mu2 * strainIncrement(2, 2) + lam * (strainIncrement(0, 0) + strainIncrement(1, 1));

	stressIncrement(0, 1) = mu * strainIncrement(0, 1);
	stressIncrement(1, 2) = mu * strainIncrement(1, 2);
	stressIncrement(0, 2) = mu * strainIncrement(0, 2);
	stressIncrement(1, 0) = mu * strainIncrement(1, 0);
	stressIncrement(2, 1) = mu * strainIncrement(2, 1);
	stressIncrement(2, 0) = mu * strainIncrement(2, 0);
}

bool MODEL::isEndingPoint() {
	double p, q, epsilonv;
	if (testType == 2 || testType == 5) {
		if (loopCounter >= loop * 2 + 1)
			return true;
		return false;
	}
	else {
		switch (endAndReversalType) {
		case 0:
			p = tr(stress) / 3;
			if (abs(p) >= abs(endAndReversalPoint)) {
				return true;
			}
			else
				return false;

		case 1:
			q = stress(0, 0) - stress(2, 2);
			if (abs(q) >= abs(endAndReversalPoint)) {
				return true;
			}
			else
				return false;

		case 2:
			epsilonv = tr(strain);
			if (abs(epsilonv) >= abs(endAndReversalPoint))
				return true;
			return false;

		default:
			return true;
		}
	}
}

bool MODEL::isReversalPoint() {
	double q;
	if (testType != 2 && testType != 5)
		return false;
	switch (endAndReversalType) {
	case 0:
		return false;

	case 1:
		q = stress(0, 0) - stress(3, 3);
		if (q >= abs(endAndReversalPoint)) {
			direction = false;
			loopCounter += 1;
			return true;
		}
		else if (q < -abs(endAndReversalPoint)) {
			direction = true;
			loopCounter += 1;
			return true;
		}
		else
			return false;

	case 2:
		return false;

	default:
		return false;
	}
}