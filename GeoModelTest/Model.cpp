#include "Model.h"

MODEL::MODEL() {
	testType = 0;
	model = 0;
	pAtmos = 101.4;

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
		strainIncrement(0, 0) = -stepLength / 2;
		strainIncrement(1, 1) = -stepLength / 2;
		strainIncrement(2, 2) = stepLength;
		break;

	case 1:
		strainIncrement.clear();
		strainIncrement(0, 0) = stepLength / 2;
		strainIncrement(1, 1) = stepLength / 2;
		strainIncrement(2, 2) = -stepLength;
		break;

	case 2:
		strainIncrement.clear();
		if (direction) {
			strainIncrement(0, 0) = -stepLength / 2;
			strainIncrement(1, 1) = -stepLength / 2;
			strainIncrement(2, 2) = stepLength;
		}
		else {
			strainIncrement(0, 0) = stepLength / 2;
			strainIncrement(1, 1) = stepLength / 2;
			strainIncrement(2, 2) = -stepLength;
		}
		break;

	case 3:
		strainIncrement.clear();
		strainIncrement(0, 0) = 0;
		strainIncrement(1, 1) = 0;
		strainIncrement(2, 2) = stepLength;
		do {
			Integrator(false);
			if (stressIncrement(0, 0) < 0) 
				strainIncrement(0, 0) -= stepLength / 100;
			else
				strainIncrement(0, 0) += stepLength / 100;
			if (stressIncrement(1, 1) < 0)
				strainIncrement(1, 1) -= stepLength / 100;
			else
				strainIncrement(1, 1) += stepLength / 100;
		} while (abs(stressIncrement(0, 0)) > 0 || abs(stressIncrement(1, 1)) > 0);
		break;

	case 4:
		strainIncrement.clear();
		strainIncrement(0, 0) = 0;
		strainIncrement(1, 1) = 0;
		strainIncrement(2, 2) = -stepLength;
		do {
			Integrator(false);
			if (stressIncrement(0, 0) < 0)
				strainIncrement(0, 0) -= stepLength / 100;
			else
				strainIncrement(0, 0) += stepLength / 100;
			if (stressIncrement(1, 1) < 0)
				strainIncrement(1, 1) -= stepLength / 100;
			else
				strainIncrement(1, 1) += stepLength / 100;
		} while (abs(stressIncrement(0, 0)) > stressTolerance || abs(stressIncrement(1, 1)) > stressTolerance);
		break;

	case 5:
		strainIncrement.clear();
		strainIncrement(0, 0) = 0;
		strainIncrement(1, 1) = 0;
		if (direction)
			strainIncrement(2, 2) = stepLength;
		else
			strainIncrement(2, 2) = -stepLength;
		do {
			Integrator(false);
			if (stressIncrement(0, 0) < 0)
				strainIncrement(0, 0) -= stepLength / 100;
			else
				strainIncrement(0, 0) += stepLength / 100;
			if (stressIncrement(1, 1) < 0)
				strainIncrement(1, 1) -= stepLength / 100;
			else
				strainIncrement(1, 1) += stepLength / 100;
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
	case 2:
		if (int(internalParameter[0]) == 1)
			IntegratorDMExplicit(updateFlag);
		else if (int(internalParameter[0]) == 0)
			IntegratorDMImplicit(updateFlag);
		break;
	default:
		stressIncrement.clear();
		break;
	}
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

void MODEL::IntegratorDMExplicit(bool updateFlag) {
	double G, K, f, B = 0, C = 0, D = 0, g = 0, L = 0;
	MATRIX ds, n, alpha, z, alphaInit, dz, dAlpha;
	for (int i = 0; i < 9; i++)
		alpha.matrix[i] = saveParameter.back().at(i + 1);
	for (int i = 0; i < 9; i++)
		z.matrix[i] = saveParameter.back().at(i + 10);
	for (int i = 0; i < 9; i++)
		alphaInit.matrix[i] = saveParameter.back().at(i + 19);
	double p = tr(stress) / 3; 
	double ee = saveParameter.back().at(0);
	
	G = getG(p, ee);
	K = getK(G);
	ds = 2 * G * strainIncrement;
	f = getF(stress + ds, alpha, p);

	if (f < 0) {
		alphaInit = alpha;
	}
	else {

	}
}

void MODEL::IntegratorDMImplicit(bool updateFlag) {

}

double MODEL::relu(double x) {
	if (x > 0)
		return x;
	else
		return 0;
}

double MODEL::getG(double p, double e) {
	if (p < 0)
		p = 0;
	return internalParameter[1] * pAtmos * pow(2.97 - e, 2) / (1 + e) * pow(p / pAtmos, 0.5);
}

double MODEL::getK(double G) {
	return 2.0 / 3.0 * (1 + internalParameter[2]) / (1 - 2 * internalParameter[2]) * G;
}

double MODEL::getEc(double pc) {
	if (pc < 0)
		pc = 1e-6;
	return internalParameter[6] - internalParameter[5] * pow(pc / pAtmos, internalParameter[7]);
}

double MODEL::getF(MATRIX s, MATRIX alpha, double p) {
	return (s - p * alpha) % (s - p * alpha) - sqrt(2.0 / 3.0) * p * internalParameter[8];
}

MATRIX MODEL::getN(MATRIX r, MATRIX alpha) {
	MATRIX n;
	n = (r - alpha) / sqrt(2.0 / 3) / internalParameter[8];
	return n;
}

double MODEL::getCos3Theta(MATRIX n) {
	double cos3Theta = sqrt(6.0) * tr(n * n * n);
	if (cos3Theta > 1.0)
		return 1.0;
	else if (cos3Theta < -1.0)
		return -1.0;
	else
		return cos3Theta;
}

double MODEL::getg(double cos3Theta) {
	return 2.0 * internalParameter[4] / (1 + internalParameter[4] - (1 - internalParameter[4]) * cos3Theta);
}

MATRIX MODEL::getAlphaThetaD(MATRIX n, double p, double e, double g) {
	if (p < 0)
		p = 1e-6;
	return sqrt(2.0 / 3) * (g * internalParameter[3] * pow(2.7182818, internalParameter[13] * (e - getEc(p))) - internalParameter[8]) * n;
}

double MODEL::getB(double cos3Theta, double g) {
	return 1 + 3.0 / 2 * (1 + internalParameter[4]) / internalParameter[4] * g * cos3Theta;
}

double MODEL::getC(double g) {
	return 3.0 * sqrt(1.5) * (1 - internalParameter[4]) / internalParameter[4] * g;
}

double MODEL::getAd(MATRIX z, MATRIX n) {
	return internalParameter[12] * (1 + relu(z % n));
}

double MODEL::getD(MATRIX n, double Ad, MATRIX alpha, MATRIX alphaThetaD) {
	return Ad * ((alphaThetaD - alpha) % n);
}

MATRIX MODEL::getAlphaThetaB(MATRIX n, double p, double e, double g) {
	if (p < 0)
		p = 1e-6;
	return sqrt(2.0 / 3) * (g * internalParameter[3] * pow(2.7182818, -internalParameter[11] * (e - getEc(p))) - internalParameter[8]) * n;
}

double MODEL::getH(double e, double p, MATRIX alpha, MATRIX alphaInit, MATRIX n) {
	double b0 = internalParameter[0] * internalParameter[9] * (1 - internalParameter[10] * e) * pow(p / pAtmos, -0.5);
	double tmp = tr((alpha - alphaInit) * n);
	if (tmp < 0)
		tmp = 1e-6;
	return b0 / tmp;
}

double MODEL::getKp(MATRIX alphaThetaB, MATRIX alpha, double p, MATRIX n, double h) {
	if (p < 0)
		p = 0;
	return 2.0 / 3 * p * h * ((alphaThetaB - alpha) % n);
}

double MODEL::getL(MATRIX n, double G, MATRIX r, MATRIX de, double depsv, double Kp, double B, double C, double K, double D) {
	return (2 * G * (n % de) - (n % r) * depsv) / (Kp + 2 * G * (B - C * tr(n * n * n)) - K * D * (n % r));
}

MATRIX MODEL::getRAp(double B, double C, MATRIX n) {
	MATRIX I(1, 1, 1);
	return B * n - C * (n * n - 1.0 / 3 * I);
}

MATRIX MODEL::getdSigma(double G, MATRIX de, double K, double depsv, double L, MATRIX RAp, double D) {
	MATRIX I(1, 1, 1);
	return 2 * G * de + K * depsv * I - relu(L) * (2 * G * RAp + K * D * I);
}

MATRIX MODEL::getdAlpha(double L, double h, MATRIX alphaThetaB, MATRIX alpha) {
	return relu(L) * 2.0 / 3.0 * h * (alphaThetaB - alpha);
}

MATRIX MODEL::getdz(double depsvp, MATRIX n, MATRIX z) {
	return -internalParameter[15] * relu(-depsvp) * (internalParameter[14]* n + z);
}