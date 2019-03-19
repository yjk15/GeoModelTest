#include "Model.h"

MODEL::MODEL() {
	testType = 0;
	model = 0;
	pAtmos = 101.4;

	ee = 0;
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
	if (model == 0) {
		vector<double> tmpPara;
		tmpPara.push_back(ee);
		saveParameter.push_back(tmpPara);
	}
	else if (model == 1) {
		vector<double> tmpPara;
		tmpPara.push_back(ee);
		tmpPara.push_back(internalParameter[4]);
		saveParameter.push_back(tmpPara);
	}
	else if (model == 2) {
		MATRIX I(1, 1, 1), alpha = (stress - tr(stress) / 3 * I) / (tr(stress) / 3), alphaInit = alpha, z;
		vector<double> tmpPara;
		tmpPara.push_back(ee);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(alpha.matrix[i]);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(z.matrix[i]);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(alphaInit.matrix[i]);
		saveParameter.push_back(tmpPara);
	}
	
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
	case 1:
		IntegratorEB(updateFlag);
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
		isReversalPoint();
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
		q =  - stress(0, 0) + stress(2, 2);
		if (q >= abs(endAndReversalPoint)) {
			if (direction != false)
				loopCounter += 1;
			direction = false;
			return true;
		}
		else if (q < -abs(endAndReversalPoint)) {
			if (direction != true)
				loopCounter += 1;
			direction = true;
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

	double depsv = tr(strainIncrement) / 3;
	double ee = saveParameter.back().at(0);
	double dee = -depsv * (1 + ee);

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

	if (updateFlag) {
		vector<double> tmpPara;
		tmpPara.push_back(ee + dee);
		saveParameter.push_back(tmpPara);
	}
}

void MODEL::IntegratorDMExplicit(bool updateFlag) {
	double G, K, f, B = 0, C = 0, D = 0, g = 0, L = 0, depsv;
	MATRIX ds, n, alpha, z, alphaInit, dz, dAlpha, I(1, 1, 1);
	for (int i = 0; i < 9; i++)
		alpha.matrix[i] = saveParameter.back().at(i + 1);
	for (int i = 0; i < 9; i++)
		z.matrix[i] = saveParameter.back().at(i + 10);
	for (int i = 0; i < 9; i++)
		alphaInit.matrix[i] = saveParameter.back().at(i + 19);
	double p = tr(stress) / 3; 
	double ee = saveParameter.back().at(0);
	depsv = tr(strainIncrement) / 3;

	G = getG(p, ee);
	K = getK(G);
	ds = 2 * G * (strainIncrement - depsv * I);
	f = getF(stress - p * I + ds, alpha, p);

	if (f < 0) {
		alphaInit = alpha;
		stressIncrement = ds + depsv * K * I;
	}
	else {
		MODEL::RK4Class rk41, rk42, rk43, rk44;
		rk41 = RK4(stress, alpha, ee, z, strain, alphaInit);
		rk42 = RK4(stress + rk41.ds / 2, alpha + rk41.dAlpha / 2, ee + rk41.dee, z + rk41.dz / 2, strain, alphaInit);
		rk43 = RK4(stress + rk42.ds / 2, alpha + rk42.dAlpha / 2, ee + rk42.dee, z + rk42.dz / 2, strain, alphaInit);
		rk44 = RK4(stress + rk43.ds, alpha + rk43.dAlpha, ee + rk43.dee, z + rk43.dz, strain, alphaInit);
		stressIncrement = (rk41.ds + rk42.ds * 2 + rk43.ds * 2 + rk44.ds) / 6;
		alpha = alpha + (rk41.dAlpha + rk42.dAlpha * 2 + rk43.dAlpha * 2 + rk44.dAlpha) / 6;
		z = z + (rk41.dz + rk42.dz * 2 + rk43.dz * 2 + rk44.dz) / 6;
		ee = ee + (rk41.dee + rk42.dee * 2 + rk43.dee * 2 + rk44.dee) / 6;
	}

	f = getF(stress - (p + tr(stressIncrement)) * I + stressIncrement, alpha, p + tr(stressIncrement));

	vector<double> tmpPara;
	tmpPara.clear();
	if (updateFlag) {
		tmpPara.push_back(ee);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(alpha.matrix[i]);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(z.matrix[i]);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(alphaInit.matrix[i]);
		saveParameter.push_back(tmpPara);
	}
}

MODEL::RK4Class MODEL::RK4(MATRIX stress, MATRIX alpha, double ee, MATRIX z, MATRIX srtain, MATRIX alphaInit) {
	MATRIX I(1, 1, 1);
	double p = tr(stress) / 3;
	double G = getG(p, ee);
	double K = getK(G);
	MATRIX s = stress - p * I;
	MATRIX r = s / p;
	MATRIX n = getN(r, alpha);
	double cos3Theta = getCos3Theta(n);
	double g = getg(cos3Theta);
	MATRIX alphaThetaD = getAlphaThetaD(n, p, ee, g);
	MATRIX alphaThetaB = getAlphaThetaB(n, p, ee, g);
	double Ad = getAd(z, n);
	double B = getB(cos3Theta, g);
	double C = getC(g);
	double D = getD(n, Ad, alpha, alphaThetaD);
	double h = getH(ee, p, alpha, alphaInit, n);
	double Kp = getKp(alphaThetaB, alpha, p, n, h);
	double depsv = tr(strain) / 3;
	MATRIX de = strainIncrement - depsv * I;
	double L = getL(n, G, r, de, depsv, Kp, B, C, K, D);
	MATRIX RAp = getRAp(B, C, n);
	
	MATRIX ds = getdSigma(G, de, K, depsv, L, RAp, D);
	MATRIX dAlpha = getdAlpha(L, h, alphaThetaB, alpha);
	MATRIX dz = getdz(relu(L) * D, n, z);
	double dee = -depsv * (1 + ee);

	MODEL::RK4Class result;
	result.dAlpha = dAlpha;
	result.dz = dz;
	result.ds = ds;
	result.dee = dee;
	return result;
}

void MODEL::IntegratorDMImplicit(bool updateFlag) {
	double G, K, f, B = 0, C = 0, D = 0, g = 0, L = 0, depsv, dp = 0;
	MATRIX r, ds, n, alpha, z, alphaInit, dz, dAlpha, I(1, 1, 1);
	for (int i = 0; i < 9; i++)
		alpha.matrix[i] = saveParameter.back().at(i + 1);
	for (int i = 0; i < 9; i++)
		z.matrix[i] = saveParameter.back().at(i + 10);
	for (int i = 0; i < 9; i++)
		alphaInit.matrix[i] = saveParameter.back().at(i + 19);
	double p = tr(stress) / 3;
	double ee = saveParameter.back().at(0);
	depsv = tr(strainIncrement) / 3;
	MATRIX s = stress - p * I;

	G = getG(p, ee);
	K = getK(G);
	ds = 2 * G * (strainIncrement - depsv * I);
	f = getF(stress - p * I + ds, alpha, p);
	MATRIX alphaThetaD, alphaThetaB, RAp, tmp;
	double cos3Theta, Ad, h, Kp, dL, dee;
	dee = -depsv * (1 + ee);

	if (f < 0) {
		alphaInit = alpha;
		stressIncrement = ds + depsv * K * I;
	}
	else {
		for (int i = 0; i < 100; i++) {
			G = getG(p + dp, ee);
			K = getK(G);
			r = (s + ds) / (p + dp);
			n = getN(r, alpha + dAlpha);
			cos3Theta = getCos3Theta(n);
			g = getg(cos3Theta);
			alphaThetaD = getAlphaThetaD(n, p + dp, ee + dee, g);
			alphaThetaB = getAlphaThetaB(n, p + dp, ee + dee, g);
			Ad = getAd(z + dz, n);
			B = getB(cos3Theta, g);
			C = getC(g);
			D = getD(n, Ad, alpha + dAlpha, alphaThetaD);
			h = getH(ee + dee, p + dp, alpha + dAlpha, alphaInit, n);
			RAp = getRAp(B, C, n);
			Kp = getKp(alphaThetaB, alpha + dAlpha, p + dp, n, h);
			f = getF(s + ds, alpha + dAlpha, p + dp);
			if (abs(f) < 1e-6)
				break;

			tmp = -2 * G * RAp + D * K * alpha - 2.0 / 3 * p * h * (alphaThetaB - alpha);
			dL = f / ((tmp % (s + ds - (p + dp) * (alpha + dAlpha))) / sqrt((s + ds - (p + dp) * (alpha + dAlpha)) % (s + ds - (p + dp) * (alpha + dAlpha))) + sqrt(2.0 / 3) * D * K * internalParameter[8]);

			L = L - dL;
			dp = K * (-relu(L) * D);
			ds = (strainIncrement - depsv * I - relu(L) * RAp) * 2 * G;
			dAlpha = getdAlpha(L, h, alphaThetaB, alpha + dAlpha);
			dz = getdz(relu(L) * D, n, z + dz);
		}
		alpha = alpha + dAlpha;
		z = z + dz;
		stressIncrement = ds + dp * I;
		ee = ee + dee;
	}

	vector<double> tmpPara;
	tmpPara.clear();
	if (updateFlag) {
		tmpPara.push_back(ee);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(alpha.matrix[i]);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(z.matrix[i]);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(alphaInit.matrix[i]);
		saveParameter.push_back(tmpPara);
	}
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
	return sqrt((s - p * alpha) % (s - p * alpha)) - sqrt(2.0 / 3.0) * p * internalParameter[8];
}

MATRIX MODEL::getN(MATRIX r, MATRIX alpha) {
	MATRIX n, I(-1/sqrt(6), -1 / sqrt(6), 2/sqrt(6));
	n = (r - alpha) / sqrt(2.0 / 3) / internalParameter[8];
	if ((n % n) - 1 < 1e-6) {
		if (abs(n(0, 0)) > 1e-6)
			n = I;
		else
			n = I * n(0, 0) / abs(n(0, 0));
	}
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
	return 1 + 3.0 / 2 * (1 - internalParameter[4]) / internalParameter[4] * g * cos3Theta;
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
	double b0 = internalParameter[1] * internalParameter[9] * (1 - internalParameter[10] * e) * pow(p / pAtmos, -0.5);
	double tmp = (alpha - alphaInit) % n;
	if (abs(tmp) < 1e-8)
		tmp = 1e-8;
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

void MODEL::IntegratorEB(bool updateFlag) {
	double Kd1 = internalParameter[0];
	double Kd2 = internalParameter[1];
	double nd2 = internalParameter[2];
	double nud = internalParameter[3];

	double depsv = tr(strainIncrement) / 3;
	double ee = saveParameter.back().at(0);
	double dee = -depsv * (1 + ee);
	double gammamax = saveParameter.back().at(1);
	double gamma = sqrt(2.0 / 3 / 3 * (pow(strain(0, 0) - strain(1, 1), 2) + pow(strain(1, 1) - strain(2, 2), 2) + pow(strain(2, 2) - strain(0, 0), 2) + 1.5*(pow(strain(0, 1), 2) + pow(strain(1, 2), 2) + pow(strain(0, 2), 2))));
	double p = 1.0 / 3 * (stress(0, 0) + stress(1, 1) + stress(2, 2));

	if (gamma > gammamax)
		gammamax = gamma;
	if (p / pAtmos < 0.01)
		p = 0.01 * pAtmos;

	double gammac = 2 * 0.65 * gammamax / pow(p / pAtmos, 1 - nd2);
	double G = Kd2 * pAtmos * sqrt(p / pAtmos) / (1 + Kd1 * gammac);
	double v = nud;
	double E = 2 * (1 + v)*G;

	double mu2 = E / (1.0 + v);
	double lam = v * mu2 / (1.0 - 2.0*v);
	double mu = 0.50*mu2;

	stressIncrement(0, 0) = mu2 * strainIncrement(0, 0) + lam * (strainIncrement(1, 1) + strainIncrement(2, 2));
	stressIncrement(1, 1) = mu2 * strainIncrement(1, 1) + lam * (strainIncrement(0, 0) + strainIncrement(2, 2));
	stressIncrement(2, 2) = mu2 * strainIncrement(2, 2) + lam * (strainIncrement(0, 0) + strainIncrement(1, 1));

	stressIncrement(0, 1) = mu * strainIncrement(0, 1);
	stressIncrement(1, 2) = mu * strainIncrement(1, 2);
	stressIncrement(0, 2) = mu * strainIncrement(0, 2);
	stressIncrement(1, 0) = mu * strainIncrement(1, 0);
	stressIncrement(2, 1) = mu * strainIncrement(2, 1);
	stressIncrement(2, 0) = mu * strainIncrement(2, 0);

	if (updateFlag) {
		vector<double> tmpPara;
		tmpPara.push_back(ee + dee);
		tmpPara.push_back(gammamax);
		saveParameter.push_back(tmpPara);
	}
}