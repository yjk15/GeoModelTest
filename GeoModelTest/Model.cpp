#include "Model.h"

MODEL::MODEL() {
	testType = 0;
	model = 0;
	pAtmos = 101.4;

	ee = 0;
	endAndReversalType = 0;
	endAndReversalPoint = 0;
	reverse = 0;
	reverseCounter = 0;
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

	timer = 0;
	betaTimer = 0;
	CPMTimer = 0;
	preTimer = 0;
	subSteps = 0;
	CPM = 0;
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
	else if (model == 3) {
		MATRIX I(1, 1, 1), alpha = (stress - tr(stress) / 3 * I) / (tr(stress) / 3);
		vector<double> tmpPara;
		tmpPara.push_back(ee);
		tmpPara.push_back(0.0);
		tmpPara.push_back(0.0);
		tmpPara.push_back(0.0);
		tmpPara.push_back(0.0);
		tmpPara.push_back(0.0);
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(alpha.matrix[i]);
		saveParameter.push_back(tmpPara);
	}

	double initPressure = stress(0,0);

	while (!isEndingPoint()) {
		stressPath->push_back(stress);
		strainPath->push_back(strain);
		GetStrainIncrementForSpecifiedTestType(initPressure);
		Integrator(true);
		strain = strain + strainIncrement;
		stress = stress + stressIncrement;
	}
	stressPath->push_back(stress);
	strainPath->push_back(strain);
}

void MODEL::GetStrainIncrementForSpecifiedTestType(double initPressure) {
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
			if (stressIncrement(0, 0) + stress(0, 0) < initPressure)
				strainIncrement(0, 0) += stepLength / 100;
			else
				strainIncrement(0, 0) -= stepLength / 100;
			if (stressIncrement(1, 1) + stress(1, 1) < initPressure)
				strainIncrement(1, 1) += stepLength / 100;
			else
				strainIncrement(1, 1) -= stepLength / 100;
		} while (abs(stressIncrement(0, 0) + stress(0, 0) - initPressure) > stressTolerance
			|| abs(stressIncrement(1, 1) + stress(1, 1) - initPressure) > stressTolerance);
		break;

	case 4:
		strainIncrement.clear();
		strainIncrement(0, 0) = 0;
		strainIncrement(1, 1) = 0;
		strainIncrement(2, 2) = -stepLength;
		do {
			Integrator(false);
			if (stressIncrement(0, 0) + stress(0, 0) < initPressure)
				strainIncrement(0, 0) += stepLength / 100;
			else
				strainIncrement(0, 0) -= stepLength / 100;
			if (stressIncrement(1, 1) + stress(1, 1) < initPressure)
				strainIncrement(1, 1) += stepLength / 100;
			else
				strainIncrement(1, 1) -= stepLength / 100;
		} while (abs(stressIncrement(0, 0) + stress(0, 0) - initPressure) > stressTolerance
			|| abs(stressIncrement(1, 1) + stress(1, 1) - initPressure) > stressTolerance);
		break;

	case 5:
		strainIncrement.clear();
		if (direction) {
			strainIncrement(2, 0) = stepLength;
			strainIncrement(0, 2) = stepLength;
		}
		else {
			strainIncrement(2, 0) = -stepLength;
			strainIncrement(0, 2) = -stepLength;
		}
		break;
	default:
		strainIncrement.clear();
		break;
	}
}

void MODEL::Integrator(bool updateFlag) {
	clock_t start, finish;
	start = clock();

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
		else if (int(internalParameter[0]) == 2)
			IntegratorDMCPPM(updateFlag);
		break;
	case 3:
		if (int(internalParameter[0]) == 1)
			IntegratorCycliqExplicit(updateFlag);
		else if (int(internalParameter[0]) == 0)
			IntegratorCycliq(updateFlag);
		break;
	default:
		stressIncrement.clear();
		break;
	}

	finish = clock();
	if (updateFlag) {
		timer += (double)(finish - start) / CLOCKS_PER_SEC;
	}

}

bool MODEL::isEndingPoint() {
	double p, q, epsilonv, epsilonq;
	if (testType == 2 || testType == 5) {
		isReversalPoint();
		if (reverseCounter >= reverse)
			return true;
		return false;
	}
	else {
		switch (endAndReversalType) {
		case 0:
			p = tr(stress) / 3;
			if (abs(abs(p) - abs(endAndReversalPoint)) < 1) {
				return true;
			}
			else
				return false;

		case 1:
			q = stress(0, 0) - stress(2, 2);
			if (abs(abs(q) - abs(endAndReversalPoint)) < 0.3) {
				return true;
			}
			else
				return false;

		case 2:
			epsilonv = tr(strain);
			if (abs(epsilonv) >= abs(endAndReversalPoint))
				return true;
			return false;

		case 3:
			epsilonq = (strain(2, 2) - strain(0, 0)) * 2.0 / 3;
			if (abs(abs(epsilonq) - abs(endAndReversalPoint)) < 1e-6)
				return true;
			return false;

		case 4:
			stepCounter += 1;
			if (stepCounter >= endAndReversalPoint)
				return true;
			return false;
		default:
			return true;
		}
	}
}

bool MODEL::isReversalPoint() {
	double q, epsq;
	if (testType != 2 && testType != 5)
		return false;
	switch (endAndReversalType) {
	case 0:
		return false;

	case 1:
		if (testType == 2)
			q = -stress(0, 0) + stress(2, 2);
		else if (testType == 5)
			q = stress(2, 0);
		if (q >= abs(endAndReversalPoint)) {
			if (direction != false)
				reverseCounter += 1;
			direction = false;
			return true;
		}
		else if (q <= -abs(endAndReversalPoint)) {
			if (direction != true)
				reverseCounter += 1;
			direction = true;
			return true;
		}
		else
			return false;

	case 2:
		return false;

	case 3:
		if (testType < 5) {
			epsq = 2.0 / 3 * (strain(2, 2) - strain(0, 0));
			if (epsq > abs(endAndReversalPoint)) {
				if (direction != false)
					reverseCounter += 1;
				direction = false;
				return true;
			}
			else if (epsq < -abs(endAndReversalPoint)) {
				if (direction != true)
					reverseCounter += 1;
				direction = true;
				return true;
			}

		}
		else if (testType == 5) {
			epsq = strain(2, 0);
			if (epsq > abs(endAndReversalPoint)) {
				if (direction != false)
					reverseCounter += 1;
				direction = false;
				return true;
			}
			else if (epsq < -abs(endAndReversalPoint)) {
				if (direction != true)
					reverseCounter += 1;
				direction = true;
				return true;
			}
		}
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
	MODEL::RK4Class rk41, rk42, rk43, rk44;
	for (int i = 0; i < 9; i++)
		alpha.matrix[i] = saveParameter.back().at(i + 1);
	for (int i = 0; i < 9; i++)
		z.matrix[i] = saveParameter.back().at(i + 10);
	for (int i = 0; i < 9; i++)
		alphaInit.matrix[i] = saveParameter.back().at(i + 19);
	double dp = 0, p = tr(stress) / 3;
	double ee = saveParameter.back().at(0);
	depsv = tr(strainIncrement);

	G = getG(p, ee);
	K = getK(G);
	ds = 2 * G * (strainIncrement - depsv / 3 * I);
	dp = K * depsv;
	f = getF(stress - p * I + ds, alpha, p + dp);

	if (f < 0) {
		alphaInit = alpha;
		stressIncrement = ds + depsv * K * I;
	}
	else {
		rk41 = RK4(stress, alpha, ee, z, strain, alphaInit);
		rk42 = RK4(stress + rk41.ds / 2, alpha + rk41.dAlpha / 2, ee + rk41.dee, z + rk41.dz / 2, strain, alphaInit);
		rk43 = RK4(stress + rk42.ds / 2, alpha + rk42.dAlpha / 2, ee + rk42.dee, z + rk42.dz / 2, strain, alphaInit);
		rk44 = RK4(stress + rk43.ds, alpha + rk43.dAlpha, ee + rk43.dee, z + rk43.dz, strain, alphaInit);
		stressIncrement = (rk41.ds + rk42.ds * 2 + rk43.ds * 2 + rk44.ds) / 6;
		alpha = alpha + (rk41.dAlpha + rk42.dAlpha * 2 + rk43.dAlpha * 2 + rk44.dAlpha) / 6;
		z = z + (rk41.dz + rk42.dz * 2 + rk43.dz * 2 + rk44.dz) / 6;
		ee = ee + (rk41.dee + rk42.dee * 2 + rk43.dee * 2 + rk44.dee) / 6;
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
	double depsv = tr(strainIncrement);
	MATRIX de = strainIncrement - depsv / 3 * I;
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
	depsv = tr(strainIncrement);
	MATRIX s = stress - p * I;
	MATRIX de = strainIncrement - depsv / 3 * I;

	G = getG(p, ee);
	K = getK(G);
	ds = 2 * G * (strainIncrement - depsv / 3 * I);
	dp = K * depsv;
	f = getF(s + ds, alpha, p + dp);
	MATRIX alphaThetaD, alphaThetaB, RAp, tmp;
	double cos3Theta, Ad, h, Kp, dL, dee;
	dee = -depsv * (1 + ee);

	if (f < 0) {
		alphaInit = alpha;
		stressIncrement = ds + depsv * K * I;
	}
	else {
		for (int i = 0; i < 100; i++) {
			if (updateFlag) {
				CPM += 1;
			}

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
			if (f < stepLength / 10)
				break;

			tmp = -2 * G * RAp + D * K * alpha - 2.0 / 3 * p * h * (alphaThetaB - alpha);
			dL = f / ((tmp % (s + ds - (p + dp) * (alpha + dAlpha))) / sqrt((s + ds - (p + dp) * (alpha + dAlpha)) % (s + ds - (p + dp) * (alpha + dAlpha))) + sqrt(2.0 / 3) * D * K * internalParameter[8]);

			L = L - dL;
			dp = K * (depsv - relu(L) * D);
			ds = (strainIncrement - depsv / 3 * I - relu(L) * RAp) * 2 * G;
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

void MODEL::IntegratorDMCPPM(bool updateFlag) {
	double G, K, f, B = 0, C = 0, D = 0, g = 0, L = 0, depsv, dp = 0, depsvp = 0;
	MATRIX r, ds, n, alpha, z, alphaInit, dz(0,0,0), dAlpha(0,0,0), I(1, 1, 1);
	for (int i = 0; i < 9; i++)
		alpha.matrix[i] = saveParameter.back().at(i + 1);
	for (int i = 0; i < 9; i++)
		z.matrix[i] = saveParameter.back().at(i + 10);
	for (int i = 0; i < 9; i++)
		alphaInit.matrix[i] = saveParameter.back().at(i + 19);
	double p = tr(stress) / 3;
	double ee = saveParameter.back().at(0);
	depsv = tr(strainIncrement);
	MATRIX s = stress - p * I;
	MATRIX de = strainIncrement - depsv / 3 * I;
	MODEL::RK4Class rk41, rk42, rk43, rk44;

	G = getG(p, ee);
	K = getK(G);
	ds = 2 * G * (strainIncrement - depsv / 3 * I);
	dp = K * depsv;
	f = getF(s + ds, alpha, p + dp);
	MATRIX alphaThetaD, alphaThetaB, RAp, dep(0,0,0), tmp, sN, alphaN, zN;
	double pN = p + dp;
	sN = s + ds;
	alphaN = alpha + dAlpha;
	zN = z + dz;
	double cos3Theta, Ad, h, Kp, dL = 0, dee;
	dee = -depsv * (1 + ee);
	double R[23], RLength = 0, RLength0 = 0, T[23][23], U[23];

	if (f < 0) {
		alphaInit = alpha;
		stressIncrement = ds + depsv * K * I;
	}
	else {
		for (int i = 0; i < 50; i++) {
			if (p + dp < 0.5)
				break;
			G = getG(p+dp, ee + dee);
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

			{
				R[0] = dp - K * (depsv - depsvp);
				R[1] = depsvp - L * D;
				tmp = ds - 2 * G * (de - dep);
				R[2] = tmp(0, 0);
				R[3] = tmp(1, 1);
				R[4] = tmp(0, 1);
				R[5] = tmp(0, 2);
				R[6] = tmp(1, 2);
				tmp = dep - L * RAp;
				R[7] = tmp(0, 0);
				R[8] = tmp(1, 1);
				R[9] = tmp(0, 1);
				R[10] = tmp(0, 2);
				R[11] = tmp(1, 2);
				tmp = dAlpha - 2.0 / 3 * L * h * (alphaThetaB - alpha);
				R[12] = tmp(0, 0);
				R[13] = tmp(1, 1);
				R[14] = tmp(0, 1);
				R[15] = tmp(0, 2);
				R[16] = tmp(1, 2);
				tmp = dz + L * internalParameter[15] * relu(-1 * D) * (internalParameter[14] * n + z + dz);
				R[17] = tmp(0, 0);
				R[18] = tmp(1, 1);
				R[19] = tmp(0, 1);
				R[20] = tmp(0, 2);
				R[21] = tmp(1, 2);
				R[22] = getF(s + ds, alpha + dAlpha, p + dp);
			}
			
			if (int(strain(2, 2)*1e6) == 208) {
				h = h;
			}

			RLength = 0;
			for (int j = 0; j < 23; j++)
				RLength += R[j] * R[j];
			RLength = sqrt(RLength);
			if (i == 0) 
				RLength0 = RLength;
			else if (i > 1) {
				if (RLength0 < RLength)
					break;
			}
			RLength0 = RLength;
			if (RLength < stepLength || R[22] < 0)
				break;

			for (int j = 0; j < 23; j++)
				for (int k = 0; k < 23; k++)
					T[j][k] = 0;
			{
				T[0][0] = 1; T[0][1] = K;
				T[1][1] = 1; T[1][22] = -D;
				T[2][2] = 1; T[2][7] = 2 * G;
				T[3][3] = 1; T[3][8] = 2 * G;
				T[4][4] = 1; T[4][9] = 2 * G;
				T[5][5] = 1; T[5][10] = 2 * G;
				T[6][6] = 1; T[6][11] = 2 * G;
				T[7][0] = -L * (-B*sN(0, 0)*sqrt(1.5) / (internalParameter[8] * pN*pN) - C * (-3.0*sN(0, 0)*(-alphaN(0, 0) + sN(0, 0) / pN) / pN / pN - 3.0*sN(0, 1)*(-alphaN(0, 1) + sN(0, 1) / pN) / pN / pN - 3.0*sN(0, 2)*(-alphaN(0, 2) + sN(0, 2) / pN) / pN / pN) / internalParameter[8] / internalParameter[8]);
				T[7][2] = -L * (sqrt(1.5)*B / (internalParameter[8] * pN) - 3.0*C*(-alphaN(0, 0) + sN(0, 0) / pN) / (pow(internalParameter[8], 2) * pN));
				T[7][4] = 3.0*C*L*(-alphaN(0, 1) + sN(0, 1) / pN) / (pow(internalParameter[8], 2) * pN);
				T[7][5] = 3.0*C*L*(-alphaN(0, 2) + sN(0, 2) / pN) / (pow(internalParameter[8], 2) * pN);
				T[7][7] = 1;
				T[7][12] = -L * (-sqrt(1.5)*B / internalParameter[8] - C * (3.0*alphaN(0, 0) - 3.0*sN(0, 0) / pN) / pow(internalParameter[8], 2));
				T[7][14] = C * L*(3.0*alphaN(0, 1) - 3.0*sN(0, 1) / pN) / pow(internalParameter[8], 2);
				T[7][15] = C * L*(3.0*alphaN(0, 2) - 3.0*sN(0, 2) / pN) / pow(internalParameter[8], 2);
				T[7][22] = -RAp(0, 0);

				T[8][0] = -L * (-sqrt(1.5)*B*sN(1, 1) / (internalParameter[8] * pN*pN) - C * (-3.0*sN(0, 1)*(-alphaN(0, 1) + sN(0, 1) / pN) / pN / pN - 3.0*sN(1, 1)*(-alphaN(1, 1) + sN(1, 1) / pN) / pN / pN - 3.0*sN(1, 2)*(-alphaN(1, 2) + sN(1, 2) / pN) / pN / pN) / pow(internalParameter[8], 2));
				T[8][3] = -L * (sqrt(1.5)*B / (internalParameter[8] * pN) - 3.0*C*(-alphaN(1, 1) + sN(1, 1) / pN) / (pow(internalParameter[8], 2) * pN));
				T[8][4] = 3.0*C*L*(-alphaN(0, 1) + sN(0, 1) / pN) / (pow(internalParameter[8], 2) * pN);
				T[8][6] = 3.0*C*L*(-alphaN(1, 2) + sN(1, 2) / pN) / (pow(internalParameter[8], 2) * pN);
				T[8][8] = 1;
				T[8][13] = -L * (-sqrt(1.5)*B / internalParameter[8] - C * (3.0*alphaN(1, 1) - 3.0*sN(1, 1) / pN) / pow(internalParameter[8], 2));
				T[8][14] = C * L*(3.0*alphaN(0, 1) - 3.0*sN(0, 1) / pN) / pow(internalParameter[8], 2);
				T[8][16] = C * L*(3.0*alphaN(1, 2) - 3.0*sN(1, 2) / pN) / pow(internalParameter[8], 2);
				T[8][22] = -RAp(1, 1);

				T[9][0] = -L * (-sqrt(1.5)*B*sN(0, 1) / (internalParameter[8] * pN*pN) - C * (-1.5*sN(0, 0)*(-alphaN(0, 1) + sN(0, 1) / pN) / pN / pN - sN(0, 1) * (-1.5*alphaN(0, 0) + 1.5*sN(0, 0) / pN) / pN / pN - 1.5*sN(0, 1)*(-alphaN(1, 1) + sN(1, 1) / pN) / pN / pN - 1.5*sN(0, 2)*(-alphaN(1, 2) + sN(1, 2) / pN) / pN / pN - sN(1, 1) * (-1.5*alphaN(0, 1) + 1.5*sN(0, 1) / pN) / pN / pN - sN(1, 2) * (-1.5*alphaN(0, 2) + 1.5*sN(0, 2) / pN) / pN / pN) / pow(internalParameter[8], 2));
				T[9][2] = 1.5*C*L*(-alphaN(0, 1) + sN(0, 1) / pN) / (pow(internalParameter[8], 2) * pN);
				T[9][3] = C * L*(-1.5*alphaN(0, 1) + 1.5*sN(0, 1) / pN) / (pow(internalParameter[8], 2) * pN);
				T[9][4] = -L * (sqrt(1.5)*B / (internalParameter[8] * pN) - C * ((-1.5*alphaN(0, 0) + 1.5*sN(0, 0) / pN) / pN + 1.5*(-alphaN(1, 1) + sN(1, 1) / pN) / pN) / pow(internalParameter[8], 2));
				T[9][5] = 1.5*C*L*(-alphaN(1, 2) + sN(1, 2) / pN) / (pow(internalParameter[8], 2) * pN);
				T[9][6] = C * L*(-1.5*alphaN(0, 2) + 1.5*sN(0, 2) / pN) / (pow(internalParameter[8], 2) * pN);
				T[9][9] = 1;
				T[9][12] = C * L*(1.5*alphaN(0, 1) - 1.5*sN(0, 1) / pN) / pow(internalParameter[8], 2);
				T[9][13] = C * L*(1.5*alphaN(0, 1) - 1.5*sN(0, 1) / pN) / pow(internalParameter[8], 2);
				T[9][14] = -L * (-sqrt(1.5)*B / internalParameter[8] - C * (1.5*alphaN(0, 0) + 1.5*alphaN(1, 1) - 1.5*sN(0, 0) / pN - 1.5*sN(1, 1) / pN) / pow(internalParameter[8], 2));
				T[9][15] = C * L*(1.5*alphaN(1, 2) - 1.5*sN(1, 2) / pN) / pow(internalParameter[8], 2);
				T[9][16] = C * L*(1.5*alphaN(0, 2) - 1.5*sN(0, 2) / pN) / pow(internalParameter[8], 2);
				T[9][22] = -B * (-sqrt(1.5)*alphaN(0, 1) + sqrt(1.5)*sN(0, 1) / pN) / internalParameter[8] + C * (1.5*(-alphaN(0, 0) + sN(0, 0) / pN)*(-alphaN(0, 1) + sN(0, 1) / pN) + 1.5*(-alphaN(0, 1) + sN(0, 1) / pN)*(-alphaN(1, 1) + sN(1, 1) / pN) + 1.5*(-alphaN(0, 2) + sN(0, 2) / pN)*(-alphaN(1, 2) + sN(1, 2) / pN)) / pow(internalParameter[8], 2);

				T[10][0] = -L * (-sqrt(1.5)*B*sN(0, 2) / (internalParameter[8] * pN*pN) - C * (-1.5*sN(0, 0)*(-alphaN(0, 2) + sN(0, 2) / pN) / pN / pN - 1.5*sN(0, 1)*(-alphaN(1, 2) + sN(1, 2) / pN) / pN / pN - sN(0, 2) * (-1.5*alphaN(0, 0) + 1.5*sN(0, 0) / pN) / pN / pN - 1.5*sN(0, 2)*(-alphaN(2, 2) + sN(2, 2) / pN) / pN / pN - sN(1, 2) * (-1.5*alphaN(0, 1) + 1.5*sN(0, 1) / pN) / pN / pN - sN(2, 2) * (-1.5*alphaN(0, 2) + 1.5*sN(0, 2) / pN) / pN / pN) / pow(internalParameter[8], 2));
				T[10][2] = 1.5*C*L*(-alphaN(0, 2) + sN(0, 2) / pN) / (pow(internalParameter[8], 2) * pN);
				T[10][4] = 1.5*C*L*(-alphaN(1, 2) + sN(1, 2) / pN) / (pow(internalParameter[8], 2) * pN);
				T[10][5] = -L * (sqrt(1.5)*B / (internalParameter[8] * pN) - C * ((-1.5*alphaN(0, 0) + 1.5*sN(0, 0) / pN) / pN + 1.5*(-alphaN(2, 2) + sN(2, 2) / pN) / pN) / pow(internalParameter[8], 2));
				T[10][6] = C * L*(-1.5*alphaN(0, 1) + 1.5*sN(0, 1) / pN) / (pow(internalParameter[8], 2) * pN);
				T[10][10] = 1;
				T[10][12] = C * L*(1.5*alphaN(0, 2) - 1.5*sN(0, 2) / pN) / pow(internalParameter[8], 2);
				T[10][14] = C * L*(1.5*alphaN(1, 2) - 1.5*sN(1, 2) / pN) / pow(internalParameter[8], 2);
				T[10][15] = -L * (-sqrt(1.5)*B / internalParameter[8] - C * (1.5*alphaN(0, 0) + 1.5*alphaN(2, 2) - 1.5*sN(0, 0) / pN - 1.5*sN(2, 2) / pN) / pow(internalParameter[8], 2));
				T[10][16] = C * L*(1.5*alphaN(0, 1) - 1.5*sN(0, 1) / pN) / pow(internalParameter[8], 2);
				T[10][22] = -B * (-sqrt(1.5)*alphaN(0, 2) + sqrt(1.5)*sN(0, 2) / pN) / internalParameter[8] + C * (1.5*(-alphaN(0, 0) + sN(0, 0) / pN)*(-alphaN(0, 2) + sN(0, 2) / pN) + 1.5*(-alphaN(0, 1) + sN(0, 1) / pN)*(-alphaN(1, 2) + sN(1, 2) / pN) + 1.5*(-alphaN(0, 2) + sN(0, 2) / pN)*(-alphaN(2, 2) + sN(2, 2) / pN)) / pow(internalParameter[8], 2);

				T[11][0] = -L * (-sqrt(1.5)*B*sN(1, 2) / (internalParameter[8] * pN*pN) - C * (-1.5*sN(0, 1)*(-alphaN(0, 2) + sN(0, 2) / pN) / pN / pN - sN(0, 2) * (-1.5*alphaN(0, 1) + 1.5*sN(0, 1) / pN) / pN / pN - 1.5*sN(1, 1)*(-alphaN(1, 2) + sN(1, 2) / pN) / pN / pN - sN(1, 2) * (-1.5*alphaN(1, 1) + 1.5*sN(1, 1) / pN) / pN / pN - 1.5*sN(1, 2)*(-alphaN(2, 2) + sN(2, 2) / pN) / pN / pN - sN(2, 2) * (-1.5*alphaN(1, 2) + 1.5*sN(1, 2) / pN) / pN / pN) / pow(internalParameter[8], 2));
				T[11][3] = 1.5*C*L*(-alphaN(1, 2) + sN(1, 2) / pN) / (pow(internalParameter[8], 2) * pN);
				T[11][4] = 1.5*C*L*(-alphaN(0, 2) + sN(0, 2) / pN) / (pow(internalParameter[8], 2) * pN);
				T[11][5] = C * L*(-1.5*alphaN(0, 1) + 1.5*sN(0, 1) / pN) / (pow(internalParameter[8], 2) * pN);
				T[11][6] = -L * (sqrt(1.5)*B / (internalParameter[8] * pN) - C * ((-1.5*alphaN(1, 1) + 1.5*sN(1, 1) / pN) / pN + 1.5*(-alphaN(2, 2) + sN(2, 2) / pN) / pN) / pow(internalParameter[8], 2));
				T[11][11] = 1;
				T[11][13] = C * L*(1.5*alphaN(1, 2) - 1.5*sN(1, 2) / pN) / pow(internalParameter[8], 2);
				T[11][14] = C * L*(1.5*alphaN(0, 2) - 1.5*sN(0, 2) / pN) / pow(internalParameter[8], 2);
				T[11][15] = C * L*(1.5*alphaN(0, 1) - 1.5*sN(0, 1) / pN) / pow(internalParameter[8], 2);
				T[11][16] = -L * (-sqrt(1.5)*B / internalParameter[8] - C * (1.5*alphaN(1, 1) + 1.5*alphaN(2, 2) - 1.5*sN(1, 1) / pN - 1.5*sN(2, 2) / pN) / pow(internalParameter[8], 2));
				T[11][22] = -B * (-sqrt(1.5)*alphaN(1, 2) + sqrt(1.5)*sN(1, 2) / pN) / internalParameter[8] + C * (1.5*(-alphaN(0, 1) + sN(0, 1) / pN)*(-alphaN(0, 2) + sN(0, 2) / pN) + 1.5*(-alphaN(1, 1) + sN(1, 1) / pN)*(-alphaN(1, 2) + sN(1, 2) / pN) + 1.5*(-alphaN(1, 2) + sN(1, 2) / pN)*(-alphaN(2, 2) + sN(2, 2) / pN)) / pow(internalParameter[8], 2);

				T[12][12] = 2.0 / 3 * L*h + 1;
				T[12][22] = -2.0 / 3 * h*(-alphaN(0, 0) + alphaThetaB(0, 0));

				T[13][13] = 2.0 / 3 * L*h + 1;
				T[13][22] = -2.0 / 3 * h*(-alphaN(1, 1) + alphaThetaB(1, 1));

				T[14][14] = 2.0 / 3 * L*h + 1;
				T[14][22] = -2.0 / 3 * h*(-alphaN(0, 1) + alphaThetaB(0, 1));

				T[15][15] = 2.0 / 3 * L*h + 1;
				T[15][22] = -2.0 / 3 * h*(-alphaN(0, 2) + alphaThetaB(0, 2));

				T[16][16] = 2.0 / 3 * L*h + 1;
				T[16][22] = -2.0 / 3 * h*(-alphaN(1, 2) + alphaThetaB(1, 2));

				T[17][0] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * sN(0, 0)*internalParameter[14] / (internalParameter[8] * pN*pN);
				T[17][2] = sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / (internalParameter[8] * pN);
				T[17][12] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / internalParameter[8];
				T[17][17] = relu(-D) * L*internalParameter[15] + 1;
				T[17][22] = relu(-D) * internalParameter[15] * (zN(0, 0) + internalParameter[14] * (-sqrt(1.5)*alphaN(0, 0) + sqrt(1.5)*sN(0, 0) / pN) / internalParameter[8]);

				T[18][0] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * sN(1, 1)*internalParameter[14] / (internalParameter[8] * pN*pN);
				T[18][3] = sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / (internalParameter[8] * pN);
				T[18][13] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / internalParameter[8];
				T[18][18] = relu(-D) * L*internalParameter[15] + 1;
				T[18][22] = relu(-D) * internalParameter[15] * (zN(1, 1) + internalParameter[14] * (-sqrt(1.5)*alphaN(1, 1) + sqrt(1.5)*sN(1, 1) / pN) / internalParameter[8]);

				T[19][0] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * sN(0, 1)*internalParameter[14] / (internalParameter[8] * pN*pN);
				T[19][4] = sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / (internalParameter[8] * pN);
				T[19][14] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / internalParameter[8];
				T[19][19] = relu(-D) * L*internalParameter[15] + 1;
				T[19][22] = relu(-D) * internalParameter[15] * (zN(0, 1) + internalParameter[14] * (-sqrt(1.5)*alphaN(0, 1) + sqrt(1.5)*sN(0, 1) / pN) / internalParameter[8]);

				T[20][0] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * sN(0, 2)*internalParameter[14] / (internalParameter[8] * pN*pN);
				T[20][5] = sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / (internalParameter[8] * pN);
				T[20][15] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / internalParameter[8];
				T[20][20] = relu(-D) * L*internalParameter[15] + 1;
				T[20][22] = relu(-D) * internalParameter[15] * (zN(0, 2) + internalParameter[14] * (-sqrt(1.5)*alphaN(0, 2) + sqrt(1.5)*sN(0, 2) / pN) / internalParameter[8]);

				T[21][0] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * sN(1, 2)*internalParameter[14] / (internalParameter[8] * pN*pN);
				T[21][6] = sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / (internalParameter[8] * p);
				T[21][16] = -sqrt(1.5)*relu(-D)*L*internalParameter[15] * internalParameter[14] / internalParameter[8];
				T[21][21] = relu(-D) * L*internalParameter[15] + 1;
				T[21][22] = relu(-D) * internalParameter[15] * (zN(1, 2) + internalParameter[14] * (-sqrt(1.5)*alphaN(1, 2) + sqrt(1.5)*sN(1, 2) / pN) / internalParameter[8]);

				double tmp1 = sqrt((sN - pN * alphaN) % (sN - pN * alphaN));

				T[22][0] = -sqrt(2.0 / 3)*internalParameter[8] + (-1.0*alphaN(0, 0)*(-alphaN(0, 0) * pN + sN(0, 0)) - 2.0*alphaN(0, 1)*(-alphaN(0, 1) * pN + sN(0, 1)) - 2.0*alphaN(0, 2)*(-alphaN(0, 2) * pN + sN(0, 2)) - 1.0*alphaN(1, 1)*(-alphaN(1, 1) * pN + sN(1, 1)) - 2.0*alphaN(1, 2)*(-alphaN(1, 2) * pN + sN(1, 2)) - 1.0*alphaN(2, 2)*(-alphaN(2, 2) * pN + sN(2, 2))) / tmp1;
				T[22][2] = (-1.0*alphaN(0, 0)*pN + 1.0*sN(0, 0)) / tmp1;
				T[22][3] = (-1.0*alphaN(1, 1)*pN + 1.0*sN(1, 1)) / tmp1;
				T[22][4] = (-2.0*alphaN(0, 1)*pN + 2.0*sN(0, 1)) / tmp1;
				T[22][5] = (-2.0*alphaN(0, 2)*pN + 2.0*sN(0, 2)) / tmp1;
				T[22][6] = (-2.0*alphaN(1, 2)*pN + 2.0*sN(1, 2)) / tmp1;
				T[22][12] = -1.0*pN*(-alphaN(0, 0) * pN + sN(0, 0)) / tmp1;
				T[22][13] = -1.0*pN*(-alphaN(1, 1) * pN + sN(1, 1)) / tmp1;
				T[22][14] = -2.0*pN*(-alphaN(0, 1) * pN + sN(0, 1)) / tmp1;
				T[22][15] = -2.0*pN*(-alphaN(0, 2) * pN + sN(0, 2)) / tmp1;
				T[22][16] = -2.0*pN*(-alphaN(1, 2) * pN + sN(1, 2)) / tmp1;
				tmp = -2 * G * RAp + D * K * alphaN - 2.0 / 3 * pN * h * (alphaThetaB - alphaN);
				T[22][22] = ((tmp % (sN - (p + dp) * (alpha + dAlpha))) / sqrt((s + ds - (p + dp) * (alpha + dAlpha)) % (s + ds - (p + dp) * (alpha + dAlpha))) + sqrt(2.0 / 3) * D * K * internalParameter[8]);
			}

			double cond = Guass(T, R, U);
			dp -= U[0]; pN = p + dp;
			depsvp -= U[1];
			ds(0, 0) -= U[2]; ds(1, 1) -= U[3]; ds(2, 2) += (U[2] + U[3]); ds(0, 1) -= U[4]; ds(0, 2) -= U[5]; ds(1, 2) -= U[6]; ds(1, 0) -= U[4]; ds(2, 0) -= U[5]; ds(2, 1) -= U[6]; sN = s + ds;
			dep(0, 0) -= U[7]; dep(1, 1) -= U[8]; dep(2, 2) += (U[7] + U[8]); dep(0, 1) -= U[9]; dep(0, 2) -= U[10]; dep(1, 2) -= U[11]; dep(1, 0) -= U[9]; dep(2, 0) -= U[10]; dep(2, 1) -= U[11];
			dAlpha(0, 0) -= U[12]; dAlpha(1, 1) -= U[13]; dAlpha(2, 2) += (U[12] + U[13]); dAlpha(0, 1) -= U[14]; dAlpha(0, 2) -= U[15]; dAlpha(1, 2) -= U[16]; dAlpha(1, 0) -= U[14]; dAlpha(2, 0) -= U[15]; dAlpha(2, 1) -= U[16]; alphaN = alpha + dAlpha;
			dz(0, 0) -= U[17]; dz(1, 1) -= U[18]; dz(2, 2) += (U[17] + U[18]); dz(0, 1) -= U[19]; dz(0, 2) -= U[20]; dz(1, 2) -= U[21]; dz(1, 0) -= U[19]; dz(2, 0) -= U[20]; dz(2, 1) -= U[21]; zN = z + dz;
			dL = -U[22]; L = L + dL;
		}

		if (p + dp < 0.5 || RLength0 < RLength) {
			rk41 = RK4(stress, alpha, ee, z, strain, alphaInit);
			rk42 = RK4(stress + rk41.ds / 2, alpha + rk41.dAlpha / 2, ee + rk41.dee, z + rk41.dz / 2, strain, alphaInit);
			rk43 = RK4(stress + rk42.ds / 2, alpha + rk42.dAlpha / 2, ee + rk42.dee, z + rk42.dz / 2, strain, alphaInit);
			rk44 = RK4(stress + rk43.ds, alpha + rk43.dAlpha, ee + rk43.dee, z + rk43.dz, strain, alphaInit);
			MATRIX dSigma = (rk41.ds + rk42.ds * 2 + rk43.ds * 2 + rk44.ds) / 6;
			dp = tr(dSigma) / 3;
			ds = dSigma - dp * I;
			dAlpha = (rk41.dAlpha + rk42.dAlpha * 2 + rk43.dAlpha * 2 + rk44.dAlpha) / 6;
			dz = (rk41.dz + rk42.dz * 2 + rk43.dz * 2 + rk44.dz) / 6;
			dee = (rk41.dee + rk42.dee * 2 + rk43.dee * 2 + rk44.dee) / 6;
			CPM += 1;
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

double MODEL::Guass(double a[23][23], double b[23], double x[23]) {
	double s[22];
	for (int i = 0; i < 23; i++) {
		s[i] = abs(a[i][0]);
		for (int j = 1; j < 23; j++)
			if (abs(a[i][j]) > s[i])
				s[i] = abs(a[i][j]);
	}

	for (int i = 0; i < 23; i++) {
		for (int j = 0; j < 23; j++) {
			a[i][j] /= s[i];
		}
		b[i] /= s[i];
	}

	//double aInverse[23][23];
	//LUP_solve_inverse(a, aInverse);
	double cond = 0;
	//cond = getNorm(a) * getNorm(aInverse);
	
	int i, j, k, n = 23;
	double c[23];    //存储初等行变换的系数，用于行的相减
	//消元的整个过程如下，总共n-1次消元过程。
	for (k = 0; k < n - 1; k++)
	{
		//求出第K次初等行变换的系数
		for (i = k + 1; i < n; i++)
			c[i] = a[i][k] / a[k][k];

		//第K次的消元计算
		for (i = k + 1; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				a[i][j] = a[i][j] - c[i] * a[k][j];
			}
			b[i] = b[i] - c[i] * b[k];
		}
	}

	//先计算出最后一个未知数；
	x[n-1] = b[n-1] / a[n-1][n-1];
	//求出每个未知数的值
	for (i = n-2; i >= 0; i--)
	{
		double sum = 0;
		for (j = i + 1; j < n; j++)
		{
			sum += a[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / a[i][i];
	}

	return cond;
}

//LUP分解
void MODEL::LUP_Descomposition(double A[23][23], double L[23][23], double U[23][23], int P[23])
{
	int row = 0;
	for (int i = 0; i < 23; i++)
	{
		P[i] = i;
	}
	for (int i = 0; i < 23 - 1; i++)
	{
		double p = 0.0;
		for (int j = i; j < 23; j++)
		{
			if (abs(A[j][i]) > p)
			{
				p = abs(A[j][i]);
				row = j;
			}
		}
		if (0 == p)
		{
			//cout<< "矩阵奇异，无法计算逆" <<endl;
			return;
		}

		//交换P[i]和P[row]
		int tmp = P[i];
		P[i] = P[row];
		P[row] = tmp;

		double tmp2 = 0.0;
		for (int j = 0; j < 23; j++)
		{
			//交换A[i][j]和 A[row][j]
			tmp2 = A[i][j];
			A[i][j] = A[row][j];
			A[row][j] = tmp2;
		}

		//以下同LU分解
		double u = A[i][i], l = 0.0;
		for (int j = i + 1; j < 23; j++)
		{
			l = A[j][i] / u;
			A[j][i] = l;
			for (int k = i + 1; k < 23; k++)
			{
				A[j][k] = A[j][k] - A[i][k] * l;
			}
		}

	}

	//构造L和U
	for (int i = 0; i < 23; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (i != j)
			{
				L[i][j] = A[i][j];
			}
			else
			{
				L[i][j] = 1;
			}
		}
		for (int k = i; k < 23; k++)
		{
			U[i][k] = A[i][k];

		}

	}
}
//LUP求解方程
void MODEL::LUP_Solve(double L[23][23], double U[23][23], int P[23], double b[23], double x[23])
{
	double y[23]{};

	//正向替换
	for (int i = 0; i < 23; i++)
	{
		y[i] = b[P[i]];
		for (int j = 0; j < i; j++)
		{
			y[i] = y[i] - L[i][j] * y[j];
		}
	}
	//反向替换
	for (int i = 23 - 1; i >= 0; i--)
	{
		x[i] = y[i];
		for (int j = 23 - 1; j > i; j--)
		{
			x[i] = x[i] - U[i][j] * x[j];
		}
		x[i] /= U[i][i];
	}

}

/*****************矩阵原地转置BEGI23********************/

/* 后继 */
int MODEL::getNext(int i, int m, int n)
{
	return (i%n)*m + i / n;
}

/* 前驱 */
int MODEL::getPre(int i, int m, int n)
{
	return (i%m)*n + i / m;
}

/* 处理以下标i为起点的环 */
void MODEL::movedata(double mtx[23][23], int i, int m, int n)
{
	double temp = *(*mtx+i); // 暂存
	int cur = i;    // 当前下标
	int pre = getPre(cur, m, n);
	while (pre != i)
	{
		*(*mtx + cur) = *(*mtx + pre);
		cur = pre;
		pre = getPre(cur, m, n);
	}
	*(*mtx + cur) = temp;
}

/* 转置，即循环处理所有环 */
void MODEL::transpose(double mtx[23][23], int m, int n)
{
	for (int i = 0; i < m*n; ++i)
	{
		int next = getNext(i, m, n);
		while (next > i) // 若存在后继小于i说明重复,就不进行下去了（只有不重复时进入while循环）
			next = getNext(next, m, n);
		if (next == i)  // 处理当前环
			movedata(mtx, i, m, n);
	}
}
/*****************矩阵原地转置E23D********************/

//LUP求逆(将每列b求出的各列x进行组装)
void MODEL::LUP_solve_inverse(double A[23][23], double inv_A[23][23])
{
	//创建矩阵A的副本，注意不能直接用A计算，因为LUP分解算法已将其改变
	double A_mirror[23][23]{};
	
	double inv_A_each[23]{};//矩阵逆的各列
	//double B[23][23]{};
	double b[23]{};//b阵为B阵的列矩阵分量

	for (int i = 0; i < 23; i++)
	{
		double L[23][23]{};
		double U[23][23]{};
		int P[23]{};

		//构造单位阵的每一列
		for (int i = 0; i < 23; i++)
		{
			b[i] = 0;
		}
		b[i] = 1;

		//每次都需要重新将A复制一份
		for (int k = 0; k < 23; k++)
			for (int j = 0; j < 23; j++)
				A_mirror[k][j] = A[k][j];


		LUP_Descomposition(A_mirror, L, U, P);

		LUP_Solve(L, U, P, b, inv_A_each);
		memcpy(*inv_A + i * 23, inv_A_each, 23 * sizeof(double));//将各列拼接起来
	}
	transpose(inv_A, 23, 23);//由于现在根据每列b算出的x按行存储，因此需转置
}

double MODEL::getNorm(double a[23][23]) {
	double norm = 0, tmp[23];
	for (int i = 0; i < 23; i++) {
		tmp[i] = 0;
		for (int j = 0; j < 23; j++) {
			tmp[i] += abs(a[i][j]);
		}
		if (tmp[i] > norm)
			norm = tmp[i];
	}
	return norm;
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
	MATRIX n, N(-1 / sqrt(6), -1 / sqrt(6), 2 / sqrt(6)), M(0, 0, 0);
	n = (r - alpha) / sqrt(2.0 / 3) / internalParameter[8];
	if (abs((n % n) - 1) > 1e-2 ) {
		n = strainIncrement / sqrt(strainIncrement % strainIncrement);
		if (abs(tr(n)) > 1e-2) {
			n = N;
		}
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
	return (2 * G * (n % de) - K * (n % r) * depsv) / (Kp + 2 * G * (B - C * tr(n * n * n)) - K * D * (n % r));
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
	return -internalParameter[15] * relu(-depsvp) * (internalParameter[14] * n + z);
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

void MODEL::IntegratorCycliq(bool updateFlag) {
	clock_t start, finish, preStart, preFinish;
	preStart = clock();

	double pmin = 0.5;
	//Cycliq Model中最小的p值
	double tolerance = 1e-4;

	double ein = saveParameter.back().at(0);
	double epsvir = saveParameter.back().at(1);
	double epsvre = saveParameter.back().at(2);
	double gammamono = saveParameter.back().at(3);
	double epsvc = saveParameter.back().at(4);
	double etam = saveParameter.back().at(5);
	MATRIX alpha;
	for (int i = 0; i < 9; i++) {
		alpha.matrix[i] = saveParameter.back().at(i + 6);
	}

	double epsvir_ns, epsvre_ns, gammamonos, epsvc_ns, etamplus1;
	MATRIX alpha_ns;
	double G0, kappa, h, M, dre1, dre2, rdr, eta, dir, lamdac, ksi, e0, np, nd;
	G0 = internalParameter[1];
	kappa = internalParameter[2];
	h = internalParameter[3];
	M = internalParameter[4];
	dre1 = internalParameter[5];
	dre2 = internalParameter[6];
	dir = internalParameter[7];
	eta = internalParameter[8];
	rdr = internalParameter[9];
	np = internalParameter[10];
	nd = internalParameter[11];
	lamdac = internalParameter[12];
	e0 = internalParameter[13];
	ksi = internalParameter[14];
	double depsv;
	depsv = tr(strainIncrement);
	double Mfc = M, Mdc = M;
	double sinphi = 3.0 * Mfc / (Mfc + 6.0);
	double tanphi = sinphi / sqrt(1.0 - sinphi * sinphi);
	double Mfo = 2 * sqrt(3.0) * tanphi / sqrt(3.0 + 4.0 * tanphi * tanphi);
	double epsvir_n = epsvir;
	double epsvre_n = epsvre;
	double epsvc_n = epsvc;

	MATRIX dev_strain, dev_strain_n, ddev_strain_p, dev_stress, dev_stress_n, normal, pass, rbar, rbar0, rbar1, r1, rd, stress_n = stress, strain_n = strain, strain_nplus1 = strain_n + strainIncrement, alpha_n = alpha, stress_pass, ZeroTensor, alpha_nplus1, r, r_nplus1;

	double phi = 0, phi_n = 0, trace = 0, trace_n = 0, dtracep = 0, iconv = 0.0, lambdamax = 0, lambdamin = 0, dlambda, chi = 0, sin3theta, beta0, beta1, Fb0, Fb1, Fb = 0, gtheta, intm, beta, ec, psi;
	double N, H, loadindex, rou, roubar, eta_nplus1, Dre_n, Dir_n, D;
	int isub, sub1, sub2, sub;

	//compute the deviatoric stress of last step
	double p_n = 1.0 / 3 * tr(stress_n);
	dev_stress_n = stress_n;
	for (int i = 0; i < 3; i++)
		dev_stress_n(i, i) -= (p_n);
	if ((p_n) < pmin)
	{
		p_n = pmin;
	}

	trace = tr(strain_nplus1);
	trace_n = tr(strain_n);
	dev_strain_n = strain_n;
	for (int i = 0; i < 3; i++)
		dev_strain_n(i, i) -= (1.0 / 3 * trace_n);
	dev_strain = strain_nplus1;
	for (int i = 0; i < 3; i++)
		dev_strain(i, i) -= (1.0 / 3 * trace);

	double en = ein;//void ratio

	// --------------(I)Initialize-------------------------------------
	alpha_nplus1 = alpha_n;

	double epsvir_nplus1 = epsvir_n;
	double epsvre_nplus1 = epsvre_n;
	double epsvc_nplus1 = epsvc_n + trace - trace_n;
	etamplus1 = etam;
	double lambda = 0.0;
	double p0 = p_n;
	double epsvc0 = -2 * kappa / (1 + ein)*(sqrt(p0 / pAtmos) - sqrt(pmin / pAtmos));
	r = dev_stress_n / (p_n);
	ec = e0 - lamdac * pow(p_n / pAtmos, ksi);
	psi = en - ec;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ddev_strain_p(i, j) = 0.0;

	// --------------Trial Before Substep-------------------------------------
	double p_nplus1;
	if (epsvc_nplus1 > epsvc0)
		p_nplus1 = pAtmos * pow(sqrt(p0 / pAtmos) + (1 + ein) / 2.0 / kappa * epsvc_nplus1, 2);
	else
		p_nplus1 = pmin;
	double K, G;
	K = (1 + ein) / kappa * pAtmos*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
	G = G0 * pAtmos*(pow((2.97 - ein), 2) / (1 + ein))*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
	dev_stress = dev_stress_n + 2.0 * G * (dev_strain - dev_strain_n);

	r_nplus1 = dev_stress / p_nplus1;
	alpha_ns = alpha_n;
	epsvir_ns = epsvir_n;
	epsvre_ns = epsvre_n;
	epsvc_ns = epsvc_n;
	gammamonos = gammamono;
	double eta_n;

	sub1 = (int)(sqrt(3.0 / 2) * sqrt((r_nplus1 - r) % (r_nplus1 - r)) / 0.05) + 1;
	sub2 = (int)(sqrt(2.0 / 3 * (dev_strain - dev_strain_n) % (dev_strain - dev_strain_n)) / 0.001) + 1;
	sub = sub1;
	if (sub2 > sub1)
		sub = sub2;
	if (sub > 100)
		sub = 100;

	preFinish = clock();
	if (updateFlag) {
		preTimer += (double)(preFinish - preStart) / CLOCKS_PER_SEC;
	}
	for (isub = 0; isub < sub; isub++)
	{

		// --------------(I)Initialize-------------------------------------
		alpha_nplus1 = alpha_ns;
		epsvir_nplus1 = epsvir_ns;
		epsvre_nplus1 = epsvre_ns;
		epsvc_nplus1 = epsvc_ns + (trace - trace_n) / sub;
		r = dev_stress_n / (p_n);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				ddev_strain_p(i, j) = 0.0;
		dtracep = 0.0;
		lambda = 0.0;
		p0 = p_n;
		epsvc0 = -2 * kappa / (1 + ein)*(sqrt(p0 / pAtmos) - sqrt(pmin / pAtmos));
		if (epsvc_nplus1 < epsvc0)
			p_nplus1 = pmin;
		else
			p_nplus1 = pAtmos * pow(sqrt(p0 / pAtmos) + (1 + ein) / 2.0 / kappa * epsvc_nplus1, 2);

		K = (1 + ein) / kappa * pAtmos*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
		G = G0 * pAtmos*(pow((2.97 - ein), 2) / (1 + ein))*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
		dev_stress = dev_stress_n + 2.0 * G*(dev_strain - dev_strain_n) / sub;
		r_nplus1 = dev_stress / p_nplus1;
		eta_n = sqrt(1.5) * sqrt((r % r));

		//		r1=r/doublecontraction(r,r);
		if ((r % r) < tolerance) {
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					r1(i, j) = 0.0;
			r1(0, 0) = -1.0 / sqrt(6.0);
			r1(1, 1) = -1.0 / sqrt(6.0);
			r1(2, 2) = 2.0 / sqrt(6.0);
		}
		else {
			r1 = r / sqrt((r % r));
		}
		pass = r1 * r1*r1;
		sin3theta = -sqrt(6.0)*tr(pass);
		if (sin3theta > 1.0)
			sin3theta = 1.0;
		else if (sin3theta < -1.0)
			sin3theta = -1.0;
		gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));

		if (eta_n / gtheta > etamplus1) {
			etamplus1 = eta_n / gtheta;
		}
		if (eta_n / gtheta > Mfc*exp(-np * psi) - tolerance) {
			etamplus1 = eta_n / gtheta;
		}

		beta = getBeta(alpha_ns, r, r1, Mfc, Mfo, np, psi, etamplus1, sin3theta);
		rbar = alpha_ns + beta * (r - alpha_ns);
		normal = rbar / sqrt(rbar % rbar);
		//normal = sqrt(1.5) * normal;
		N = r % normal;

		// --------------(III)Loading/Unloading-------------------------------------

		phi = ((dev_stress - dev_stress_n) % normal) - (p_nplus1 - p_n)*N;
		phi_n = (r_nplus1 - r) % normal;
		// --------------(IV)Unloading-------------------------------------
		if (phi < tolerance || phi_n < tolerance) {
			gammamonos = 0.0;
			alpha_nplus1 = r;
			//epsvirpr=epsvir_n;
		}
		// --------------(V)Loading-------------------------------------
		else if (phi > tolerance&&phi_n > tolerance) {
			epsvc_nplus1 = epsvc_ns + (trace - trace_n) / sub - dtracep;
			loadindex = 0.0;
			lambda = 0.0;
			lambdamin = 0.0;
			lambdamax = 0.0;
			rou = sqrt(1.5) * sqrt((r - alpha_ns) % (r - alpha_ns));

			roubar = sqrt(1.5) * sqrt((rbar - alpha_ns) % (rbar - alpha_ns));

			if (roubar > tolerance) {
				H = 2.0 / 3 * h*G*gtheta*exp(-np * psi)*(Mfc*exp(-np * psi) / (etamplus1 + tolerance) *roubar / (rou + tolerance) - 1.0);
				if (H < tolerance && H >= 0) {
					H = tolerance;
				}
				if (H > -tolerance && H < 0) {
					H = -tolerance;
				}
				eta_nplus1 = sqrt(1.5) * sqrt(r_nplus1 % r_nplus1);
				rd = Mdc * exp(nd*psi) / etamplus1 * rbar;
				Dre_n = dre1 * sqrt(2.0 / 3)*((rd - r) % normal);
				if (epsvir_ns > tolerance)
					chi = -dir * epsvre_ns / epsvir_ns;
				else
					chi = 0.0;
				if (chi > 1.)
					chi = 1.;
				if (Dre_n > 0.0) {
					Dre_n = pow(-dre2 * chi, 2) / p_n;
					if (-epsvre_ns < tolerance)
						Dre_n = 0.0;
				}
				if (Dre_n > 0) {
					if (psi >= 0)
						Dir_n = dir * exp(nd*psi - eta * epsvir_ns)*(sqrt(2.0 / 3)*((rd - r) % normal))*exp(chi);
					else
						Dir_n = dir * exp(nd*psi - eta * epsvir_ns)*(sqrt(2.0 / 3)*((rd - r) % normal)*exp(chi) + pow(rdr*(1 - exp(nd*psi)) / (rdr*(1 - exp(nd*psi)) + gammamonos), 2));
				}
				else {
					if (psi >= 0)
						Dir_n = 0.0;
					else
						Dir_n = dir * exp(nd*psi - eta * epsvir_ns)*(pow(rdr*(1 - exp(nd*psi)) / (rdr*(1 - exp(nd*psi)) + gammamonos), 2));
				}

				D = Dir_n + Dre_n;

				if (strain(2, 0) > 0.0005) {
					D = D;
				}

				int wr = 1;
				iconv = 0.0;

				start = clock(); //compute the timer in the cycle
				
				do
				{
					if (iconv == 0.0) {
						//dlambda=phi/(H+2*G-K*D*N);
						dlambda = phi / abs(H + 2 * G - K * D*N);
						lambda += dlambda;
					}
					else {
						lambda = 0.5*(lambdamax + lambdamin);
					}
					loadindex = H * lambda;
					dtracep = lambda * D;
					ddev_strain_p = lambda * normal;
					epsvc_nplus1 = epsvc_ns + (trace - trace_n) / sub - dtracep;
					if (epsvc_nplus1 < epsvc0) {
						p_nplus1 = pmin;
						epsvc_nplus1 = epsvc0;
					}
					else {
						p_nplus1 = pAtmos * pow(sqrt(p0 / pAtmos) + (1 + ein) / 2.0 / kappa * epsvc_nplus1, 2);

					}
					G = G0 * pAtmos*(pow((2.97 - ein), 2) / (1 + ein))*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
					dev_stress = dev_stress_n + 2 * G*((dev_strain - dev_strain_n) / sub - ddev_strain_p);

					phi = ((dev_stress - dev_stress_n) % normal) - (p_nplus1 - p_n)*N - loadindex;
					//phi=doublecontraction(dev_stress-dev_stress_n,normal) / sqrt(1.5) -(p_nplus1-p_n)*N / sqrt(1.5)-loadindex;
					if (phi < -tolerance) {
						iconv = 1.0;
						lambdamax = lambda;
					}
					if (phi > tolerance && iconv == 1.0)
						lambdamin = lambda;
					wr = wr + 1;
					epsvir_nplus1 = lambda * Dir_n + epsvir_ns;
					epsvre_nplus1 = lambda * Dre_n + epsvre_ns;

				} while (abs(phi) > tolerance);

				finish = clock();
				CPMTimer += (double)(finish - start) / CLOCKS_PER_SEC;
				CPM += wr;
				//cout << H << "\t" << G << "\t" << K << "\t" << D << "\t" << N << endl;
				gammamonos = gammamonos + lambda;
			}

		}

		r_nplus1 = dev_stress / p_nplus1;
		eta_nplus1 = sqrt(1.5) * sqrt(r_nplus1 % r_nplus1);
		if (eta_nplus1 >= Mfc * exp(-np * psi) / (1.0 + Mfc / 3.0) - tolerance) {
			r1 = r_nplus1 / sqrt(r_nplus1 % r_nplus1);
			pass = r1 * r1*r1;
			sin3theta = -sqrt(6.0)*tr(pass);
			if (sin3theta > 1.0)
				sin3theta = 1.0;
			else if (sin3theta < -1.0)
				sin3theta = -1.0;
			gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));
			r1 = sqrt(2.0 / 3) * Mfc*exp(-np * psi)*gtheta*r1;
			if ((r1 % r1) - (r_nplus1 % r_nplus1) < tolerance) {
				intm = sqrt((r_nplus1 % r_nplus1)) / sqrt(r1 % r1) + tolerance;
				dev_stress = dev_stress / intm;
			}
			r_nplus1 = dev_stress / p_nplus1;
			eta_nplus1 = sqrt(1.5) * sqrt(r_nplus1 % r_nplus1);
		}

		alpha_ns = alpha_nplus1;
		epsvir_ns = epsvir_nplus1;
		epsvre_ns = epsvre_nplus1;

		epsvc_ns = epsvc_ns - epsvc_nplus1 + (trace - trace_n) / sub - dtracep;
		epsvc_nplus1 = epsvc_ns;
		p_n = p_nplus1;
		dev_stress_n = dev_stress;
	}
	stress_pass = dev_stress;
	for (int i = 0; i < 3; i++) {
		stress_pass(i, i) += (p_n);
	}

	stressIncrement = stress_pass - stress; // from positive to negative

	subSteps += sub;

	ein -= depsv * (1 + en);
	vector<double> tmpPara;
	if (updateFlag) {
		tmpPara.push_back(ein);
		tmpPara.push_back(epsvir_ns);
		tmpPara.push_back(epsvre_ns);
		tmpPara.push_back(gammamonos);
		tmpPara.push_back(epsvc_ns);
		tmpPara.push_back(etamplus1);
		alpha = alpha_ns;
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(alpha.matrix[i]);
		saveParameter.push_back(tmpPara);
	}
}

void MODEL::IntegratorCycliqExplicit(bool updateFlag) {
	clock_t preStart, preFinish;
	preStart = clock();

	double pmin = 0.5;
	//Cycliq Model中最小的p值
	double tolerance = 1e-4;

	double ein = saveParameter.back().at(0);
	double epsvir = saveParameter.back().at(1);
	double epsvre = saveParameter.back().at(2);
	double gammamono = saveParameter.back().at(3);
	double epsvc = saveParameter.back().at(4);
	double etam = saveParameter.back().at(5);
	MATRIX alpha;
	for (int i = 0; i < 9; i++) {
		alpha.matrix[i] = saveParameter.back().at(i + 6);
	}

	double epsvir_ns, epsvre_ns, gammamonos, epsvc_ns, etamplus1;
	MATRIX alpha_ns, I;
	double G0, kappa, h, M, dre1, dre2, rdr, eta, dir, lamdac, ksi, e0, np, nd;
	G0 = internalParameter[1];
	kappa = internalParameter[2];
	h = internalParameter[3];
	M = internalParameter[4];
	dre1 = internalParameter[5];
	dre2 = internalParameter[6];
	dir = internalParameter[7];
	eta = internalParameter[8];
	rdr = internalParameter[9];
	np = internalParameter[10];
	nd = internalParameter[11];
	lamdac = internalParameter[12];
	e0 = internalParameter[13];
	ksi = internalParameter[14];
	I(0, 0) = 1.0;
	I(1, 1) = 1.0;
	I(2, 2) = 1.0;
	double depsv;
	depsv = tr(strainIncrement);
	double Mfc = M, Mdc = M;
	double sinphi = 3.0 * Mfc / (Mfc + 6.0);
	double tanphi = sinphi / sqrt(1.0 - sinphi * sinphi);
	double Mfo = 2 * sqrt(3.0) * tanphi / sqrt(3.0 + 4.0 * tanphi * tanphi);
	double epsvir_n = epsvir;
	double epsvre_n = epsvre;
	double epsvc_n = epsvc;

	MATRIX dev_strain, dev_strain_n, ddev_strain_p, dev_stress, dev_stress_n, normal, pass, rbar, rbar0, rbar1, r1, rd, stress_n = stress, strain_n = strain, strain_nplus1 = strain_n + strainIncrement, alpha_n = alpha, stress_pass, ZeroTensor, alpha_nplus1, r, r_nplus1;

	double phi = 0, phi_n = 0, trace = 0, trace_n = 0, dtracep = 0, iconv = 0.0, lambdamax = 0, lambdamin = 0, dlambda, chi = 0, sin3theta, beta0, beta1, Fb0, Fb1, Fb = 0, gtheta, intm, beta, ec, psi;
	double N, H, loadindex, rou, roubar, eta_nplus1, Dre_n, Dir_n, D;
	int isub, sub1, sub2, sub;

	//compute the deviatoric stress of last step
	double p_n = 1.0 / 3 * tr(stress_n);
	dev_stress_n = stress_n;
	for (int i = 0; i < 3; i++)
		dev_stress_n(i, i) -= (p_n);
	if ((p_n) < pmin)
		p_n = pmin;

	trace = tr(strain_nplus1);
	trace_n = tr(strain_n);
	dev_strain_n = strain_n;
	for (int i = 0; i < 3; i++)
		dev_strain_n(i, i) -= (1.0 / 3 * trace_n);
	dev_strain = strain_nplus1;
	for (int i = 0; i < 3; i++)
		dev_strain(i, i) -= (1.0 / 3 * trace);

	double en = ein;//void ratio

	// --------------(I)Initialize-------------------------------------
	alpha_nplus1 = alpha_n;

	double epsvir_nplus1 = epsvir_n;
	double epsvre_nplus1 = epsvre_n;
	double epsvc_nplus1 = epsvc_n + trace - trace_n;
	etamplus1 = etam;
	double lambda = 0.0;
	double p0 = p_n;
	double epsvc0 = -2 * kappa / (1 + ein)*(sqrt(p0 / pAtmos) - sqrt(pmin / pAtmos));
	r = dev_stress_n / (p_n);
	ec = e0 - lamdac * pow(p_n / pAtmos, ksi);
	psi = en - ec;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ddev_strain_p(i, j) = 0.0;

	// --------------Trial Before Substep-------------------------------------
	double p_nplus1;
	if (epsvc_nplus1 > epsvc0)
		p_nplus1 = pAtmos * pow(sqrt(p0 / pAtmos) + (1 + ein) / 2.0 / kappa * epsvc_nplus1, 2);
	else
		p_nplus1 = pmin;
	double K, G;
	K = (1 + ein) / kappa * pAtmos*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
	G = G0 * pAtmos*(pow((2.97 - ein), 2) / (1 + ein))*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
	dev_stress = dev_stress_n + 2.0 * G * (dev_strain - dev_strain_n);

	r_nplus1 = dev_stress / p_nplus1;
	sub = 1;
	alpha_ns = alpha_n;
	epsvir_ns = epsvir_n;
	epsvre_ns = epsvre_n;
	epsvc_ns = epsvc_n;
	gammamonos = gammamono;
	double eta_n;

	preFinish = clock();
	if (updateFlag) {
		preTimer += (double)(preFinish - preStart) / CLOCKS_PER_SEC;
	}

	// --------------(I)Initialize-------------------------------------
	alpha_nplus1 = alpha_ns;
	epsvir_nplus1 = epsvir_ns;
	epsvre_nplus1 = epsvre_ns;
	epsvc_nplus1 = epsvc_ns + (trace - trace_n) / sub;
	r = dev_stress_n / (p_n);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ddev_strain_p(i, j) = 0.0;
	dtracep = 0.0;
	lambda = 0.0;
	p0 = p_n;
	epsvc0 = -2 * kappa / (1 + ein)*(sqrt(p0 / pAtmos) - sqrt(pmin / pAtmos));
	if (epsvc_nplus1 < epsvc0)
		p_nplus1 = pmin;
	else
		p_nplus1 = pAtmos * pow(sqrt(p0 / pAtmos) + (1 + ein) / 2.0 / kappa * epsvc_nplus1, 2);

	K = (1 + ein) / kappa * pAtmos*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
	G = G0 * pAtmos*(pow((2.97 - ein), 2) / (1 + ein))*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
	dev_stress = dev_stress_n + 2.0 * G*(dev_strain - dev_strain_n) / sub;
	r_nplus1 = dev_stress / p_nplus1;
	eta_n = sqrt(1.5) * sqrt((r % r));

	//		r1=r/doublecontraction(r,r);
	if ((r % r) < tolerance) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				r1(i, j) = 0.0;
		r1(0, 0) = -1.0 / sqrt(6.0);
		r1(1, 1) = -1.0 / sqrt(6.0);
		r1(2, 2) = 2.0 / sqrt(6.0);
	}
	else
		r1 = r / sqrt((r % r));

	pass = r1 * r1*r1;
	sin3theta = -sqrt(6.0)*tr(pass);
	if (sin3theta > 1.0)
		sin3theta = 1.0;
	else if (sin3theta < -1.0)
		sin3theta = -1.0;
	gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));

	if (eta_n / gtheta > etamplus1) {
		etamplus1 = eta_n / gtheta;
	}
	if (eta_n / gtheta > Mfc*exp(-np * psi) - tolerance) {
		etamplus1 = eta_n / gtheta;
	}

	beta = getBeta(alpha_ns, r, r1, Mfc, Mfo, np, psi, etamplus1, sin3theta);
	rbar = alpha_ns + beta * (r - alpha_ns);
	normal = rbar / sqrt(rbar % rbar);
	//normal = sqrt(1.5) * normal;
	N = r % normal;

	// --------------(III)Loading/Unloading-------------------------------------

	phi = ((dev_stress - dev_stress_n) % normal) - (p_nplus1 - p_n)*N;
	phi_n = (r_nplus1 - r) % normal;
	// --------------(IV)Unloading-------------------------------------
	if (phi < tolerance || phi_n < tolerance) {
		gammamonos = 0.0;
		alpha_ns = r;
		stressIncrement = dev_stress + p_nplus1 * I - stress;
	}
	else if (phi > tolerance && phi_n > tolerance) {
		RK4CycliqClass rk1, rk2, rk3, rk4;
		rk1 = RK4Cycliq(stress, strain, ein, epsvir, epsvre, gammamono, epsvc, etam, alpha);
		rk2 = RK4Cycliq(stress + rk1.dSigma / 2, strain, ein + rk1.dein / 2, epsvir + rk1.depsvre / 2, epsvre + rk1.depsvre / 2, gammamono + rk1.dgammamono / 2, epsvc + rk1.depsvc / 2, etam + rk1.deta / 2, alpha);
		rk3 = RK4Cycliq(stress + rk2.dSigma / 2, strain, ein + rk2.dein / 2, epsvir + rk2.depsvre / 2, epsvre + rk2.depsvre / 2, gammamono + rk2.dgammamono / 2, epsvc + rk2.depsvc / 2, etam + rk2.deta / 2, alpha);
		rk4 = RK4Cycliq(stress + rk3.dSigma, strain, ein + rk3.dein, epsvir + rk3.depsvre, epsvre + rk3.depsvre, gammamono + rk3.dgammamono, epsvc + rk3.depsvc, etam + rk3.deta, alpha);
		ein += (rk1.dein + rk2.dein * 2 + rk3.dein * 2 + rk4.dein) / 6;
		epsvir_ns = epsvir + (rk1.depsvir + rk2.depsvir * 2 + rk3.depsvir * 2 + rk4.depsvir) / 6;
		epsvre_ns = epsvre + (rk1.depsvre + rk2.depsvre * 2 + rk3.depsvre * 2 + rk4.depsvre) / 6;
		epsvc_ns = epsvc + (rk1.depsvc + rk2.depsvc * 2 + rk3.depsvc * 2 + rk4.depsvc) / 6;
		gammamonos = gammamono + (rk1.dgammamono + rk2.dgammamono * 2 + rk3.dgammamono * 2 + rk4.dgammamono) / 6;
		etamplus1 = etam + (rk1.deta + rk2.deta * 2 + rk3.deta * 2 + rk4.deta) / 6;
		stressIncrement = (rk1.dSigma + rk2.dSigma * 2 + rk3.dSigma * 2 + rk4.dSigma) / 6;
	}

	vector<double> tmpPara;
	if (updateFlag) {
		tmpPara.push_back(ein);
		tmpPara.push_back(epsvir_ns);
		tmpPara.push_back(epsvre_ns);
		tmpPara.push_back(gammamonos);
		tmpPara.push_back(epsvc_ns);
		tmpPara.push_back(etamplus1);
		alpha = alpha_ns;
		for (int i = 0; i < 9; i++)
			tmpPara.push_back(alpha.matrix[i]);
		saveParameter.push_back(tmpPara);
	}
}

MODEL::RK4CycliqClass MODEL::RK4Cycliq(MATRIX stress, MATRIX strain, double ein, double epsvir, double epsvre, double gammamono, double epsvc, double etam, MATRIX alpha) {
	clock_t preStart, preFinish;
	preStart = clock();
	
	
	double pmin = 0.5;
	//Cycliq Model中最小的p值
	double tolerance = 1e-4;

	double epsvir_ns, epsvre_ns, gammamonos, epsvc_ns, etamplus1;
	MATRIX alpha_ns, I;
	double G0, kappa, h, M, dre1, dre2, rdr, eta, dir, lamdac, ksi, e0, np, nd;
	G0 = internalParameter[1];
	kappa = internalParameter[2];
	h = internalParameter[3];
	M = internalParameter[4];
	dre1 = internalParameter[5];
	dre2 = internalParameter[6];
	dir = internalParameter[7];
	eta = internalParameter[8];
	rdr = internalParameter[9];
	np = internalParameter[10];
	nd = internalParameter[11];
	lamdac = internalParameter[12];
	e0 = internalParameter[13];
	ksi = internalParameter[14];
	double depsv;
	depsv = tr(strainIncrement);
	double Mfc = M, Mdc = M;
	double sinphi = 3.0 * Mfc / (Mfc + 6.0);
	double tanphi = sinphi / sqrt(1.0 - sinphi * sinphi);
	double Mfo = 2 * sqrt(3.0) * tanphi / sqrt(3.0 + 4.0 * tanphi * tanphi);
	double epsvir_n = epsvir;
	double epsvre_n = epsvre;
	double epsvc_n = epsvc;
	I(0, 0) = 1;
	I(1, 1) = 1;
	I(2, 2) = 1;

	if (epsvc < -1e-29) {
		I = I;
	}

	MATRIX dev_strain, dev_strain_n, ddev_strain_p, dev_stress, dev_stress_n, normal, pass, rbar, rbar0, rbar1, r1, rd, stress_n = stress, strain_n = strain, strain_nplus1 = strain_n + strainIncrement, alpha_n = alpha, stress_pass, ZeroTensor, alpha_nplus1, r, r_nplus1;

	double phi = 0, phi_n = 0, trace = 0, trace_n = 0, dtracep = 0, iconv = 0.0, lambdamax = 0, lambdamin = 0, dlambda, chi = 0, sin3theta, beta0, beta1, Fb0, Fb1, Fb = 0, gtheta, intm, beta, ec, psi;
	double N, H, loadindex, rou, roubar, eta_nplus1, Dre_n, Dir_n, D;
	int isub, sub1, sub2, sub;

	//compute the deviatoric stress of last step
	double p_n = 1.0 / 3 * tr(stress_n);
	dev_stress_n = stress_n;
	for (int i = 0; i < 3; i++)
		dev_stress_n(i, i) -= (p_n);
	if ((p_n) < pmin)
	{
		p_n = pmin;
	}

	trace = tr(strain_nplus1);
	trace_n = tr(strain_n);
	dev_strain_n = strain_n;
	for (int i = 0; i < 3; i++)
		dev_strain_n(i, i) -= (1.0 / 3 * trace_n);
	dev_strain = strain_nplus1;
	for (int i = 0; i < 3; i++)
		dev_strain(i, i) -= (1.0 / 3 * trace);

	double en = ein;//void ratio

	// --------------(I)Initialize-------------------------------------
	alpha_nplus1 = alpha_n;

	double epsvir_nplus1 = epsvir_n;
	double epsvre_nplus1 = epsvre_n;
	double epsvc_nplus1 = epsvc_n + trace - trace_n;
	etamplus1 = etam;
	double lambda = 0.0;
	double p0 = p_n;
	double epsvc0 = -2 * kappa / (1 + ein)*(sqrt(p0 / pAtmos) - sqrt(pmin / pAtmos));
	r = dev_stress_n / (p_n);
	ec = e0 - lamdac * pow(p_n / pAtmos, ksi);
	psi = en - ec;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ddev_strain_p(i, j) = 0.0;

	// --------------Trial Before Substep-------------------------------------
	double p_nplus1 = p_n;
	double K, G;
	K = (1 + ein) / kappa * pAtmos*sqrt(p_n / pAtmos);
	G = G0 * pAtmos*(pow((2.97 - ein), 2) / (1 + ein))*sqrt(p_n  / pAtmos);
	dev_stress = dev_stress_n + 2.0 * G * (dev_strain - dev_strain_n);

	r_nplus1 = dev_stress / p_nplus1;
	alpha_ns = alpha_n;
	epsvir_ns = epsvir_n;
	epsvre_ns = epsvre_n;
	epsvc_ns = epsvc_n;
	gammamonos = gammamono;
	double eta_n;

	sub1 = (int)(sqrt(3.0 / 2) * sqrt((r_nplus1 - r) % (r_nplus1 - r)) / 0.05) + 1;
	sub2 = (int)(sqrt(2.0 / 3 * (dev_strain - dev_strain_n) % (dev_strain - dev_strain_n)) / 0.001) + 1;
	sub = sub1;
	if (sub2 > sub1)
		sub = sub2;
	if (sub > 100)
		sub = 100;
	preFinish = clock();
	preTimer += (double)(preFinish - preStart) / CLOCKS_PER_SEC;
	for (isub = 0; isub < sub; isub++)
	{

		// --------------(I)Initialize-------------------------------------
		alpha_nplus1 = alpha_ns;
		epsvir_nplus1 = epsvir_ns;
		epsvre_nplus1 = epsvre_ns;
		epsvc_nplus1 = epsvc_ns + (trace - trace_n) / sub;
		r = dev_stress_n / (p_n);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				ddev_strain_p(i, j) = 0.0;
		dtracep = 0.0;
		lambda = 0.0;
		p0 = p_n;
		epsvc0 = -2 * kappa / (1 + ein)*(sqrt(p0 / pAtmos) - sqrt(pmin / pAtmos));
		dev_stress = dev_stress_n + 2.0 * G*(dev_strain - dev_strain_n) / sub;
		r_nplus1 = dev_stress / p_nplus1;
		eta_n = sqrt(1.5) * sqrt((r % r));

		//		r1=r/doublecontraction(r,r);
		if ((r % r) < tolerance) {
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					r1(i, j) = 0.0;
			r1(0, 0) = -1.0 / sqrt(6.0);
			r1(1, 1) = -1.0 / sqrt(6.0);
			r1(2, 2) = 2.0 / sqrt(6.0);
		}
		else {
			r1 = r / sqrt((r % r));
		}
		pass = r1 * r1*r1;
		sin3theta = -sqrt(6.0)*tr(pass);
		if (sin3theta > 1.0)
			sin3theta = 1.0;
		else if (sin3theta < -1.0)
			sin3theta = -1.0;
		gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));

		if (eta_n / gtheta > etamplus1)
			etamplus1 = eta_n / gtheta;
		if (eta_n / gtheta > Mfc*exp(-np * psi) - tolerance)
			etamplus1 = eta_n / gtheta;

		beta = getBeta(alpha_ns, r, r1, Mfc, Mfo, np, psi, etamplus1, sin3theta);
		rbar = alpha_ns + beta * (r - alpha_ns);
		normal = rbar / sqrt(rbar % rbar);
		//normal = sqrt(1.5) * normal;
		N = r % normal;

		rou = sqrt(1.5) * sqrt((r - alpha_ns) % (r - alpha_ns));
		roubar = sqrt(1.5) * sqrt((rbar - alpha_ns) % (rbar - alpha_ns));
		H = 2.0 / 3 * h*G*gtheta*exp(-np * psi)*(Mfc*exp(-np * psi) / (etamplus1 + tolerance) *roubar / (rou + tolerance) - 1.0);
		if (H < tolerance && H >= 0) {
			H = tolerance;
		}
		if (H > -tolerance && H < 0) {
			H = -tolerance;
		}
		eta_nplus1 = sqrt(1.5) * sqrt(r_nplus1 % r_nplus1);
		rd = Mdc * exp(nd*psi) / etamplus1 * rbar;
		Dre_n = dre1 * sqrt(2.0 / 3)*((rd - r) % normal);
		if (epsvir_ns > tolerance)
			chi = -dir * epsvre_ns / epsvir_ns;
		else
			chi = 0.0;
		if (chi > 1.)
			chi = 1.;
		if (Dre_n > 0.0) {
			Dre_n = pow(-dre2 * chi, 2) / p_n;
			if (-epsvre_ns < tolerance)
				Dre_n = 0.0;
		}
		if (Dre_n > 0) {
			if (psi >= 0) {
				Dir_n = dir * exp(nd*psi - eta * epsvir_ns)*(sqrt(2.0 / 3)*((rd - r) % normal))*exp(chi);
			}
			else {
				Dir_n = dir * exp(nd*psi - eta * epsvir_ns)*(sqrt(2.0 / 3)*((rd - r) % normal)*exp(chi) + pow(rdr*(1 - exp(nd*psi)) / (rdr*(1 - exp(nd*psi)) + gammamonos), 2));
			}
		}
		else {
			if (psi >= 0) {
				Dir_n = 0.0;
			}
			else {
				Dir_n = dir * exp(nd*psi - eta * epsvir_ns)*(pow(rdr*(1 - exp(nd*psi)) / (rdr*(1 - exp(nd*psi)) + gammamonos), 2));
			}
		}

		D = Dir_n + Dre_n;

		lambda = (2 * G * ((strainIncrement - depsv / 3 * I) % normal) - K * N * depsv) / sub / (H + 2 * G - K * D * N);
		dtracep = lambda * D;
		ddev_strain_p = lambda * normal;
		epsvc_nplus1 = epsvc_ns + (trace - trace_n) / sub - dtracep;
		if (epsvc_nplus1 < epsvc0) {
			p_nplus1 = pmin;
			epsvc_nplus1 = epsvc0;
		}
		else {
			p_nplus1 = pAtmos * pow(sqrt(p0 / pAtmos) + (1 + ein) / 2.0 / kappa * epsvc_nplus1, 2);
		}
		G = G0 * pAtmos*(pow((2.97 - ein), 2) / (1 + ein))*sqrt((p_n + p_nplus1) / 2.0 / pAtmos);
		dev_stress = dev_stress_n + 2 * G*((dev_strain - dev_strain_n) - ddev_strain_p);
		epsvir_nplus1 = lambda * Dir_n + epsvir_ns;
		epsvre_nplus1 = lambda * Dre_n + epsvre_ns;
		gammamonos = gammamonos + lambda;

		r_nplus1 = dev_stress / p_nplus1;
		eta_nplus1 = sqrt(1.5) * sqrt(r_nplus1 % r_nplus1);
		if (eta_nplus1 >= Mfc * exp(-np * psi) / (1.0 + Mfc / 3.0) - tolerance) {
			r1 = r_nplus1 / sqrt(r_nplus1 % r_nplus1);
			pass = r1 * r1*r1;
			sin3theta = -sqrt(6.0)*tr(pass);
			if (sin3theta > 1.0)
				sin3theta = 1.0;
			else if (sin3theta < -1.0)
				sin3theta = -1.0;
			gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));
			r1 = sqrt(2.0 / 3) * Mfc*exp(-np * psi)*gtheta*r1;
			if ((r1 % r1) - (r_nplus1 % r_nplus1) < tolerance) {
				intm = sqrt((r_nplus1 % r_nplus1)) / sqrt(r1 % r1) + tolerance;
				dev_stress = dev_stress / intm;
			}
			r_nplus1 = dev_stress / p_nplus1;
			eta_nplus1 = sqrt(1.5) * sqrt(r_nplus1 % r_nplus1);
		}
		alpha_ns = alpha_nplus1;
		epsvir_ns = epsvir_nplus1;
		epsvre_ns = epsvre_nplus1;
		epsvc_ns = epsvc_ns - epsvc_nplus1 + (trace - trace_n) / sub - dtracep;
		epsvc_nplus1 = epsvc_ns;
		p_n = p_nplus1;
		dev_stress_n = dev_stress;
	}

	stress_pass = dev_stress;
	for (int i = 0; i < 3; i++) {
		stress_pass(i, i) += (p_n);
	}

	subSteps += sub;

	RK4CycliqClass rk;
	rk.dSigma = stress_pass - stress;
	rk.dein = -depsv * (1 + en);
	rk.depsvc = epsvc_nplus1 - epsvc;
	rk.depsvir = epsvir_nplus1 - epsvir;
	rk.depsvre = epsvre_nplus1 - epsvre;
	rk.deta = etamplus1 - etam;
	rk.dgammamono = lambda;
	return rk;
}

double MODEL::getBeta(MATRIX alpha_ns, MATRIX r, MATRIX r1, double Mfc, double Mfo, double np, double psi, double etamplus1, double sin3theta) {
	clock_t start, finish;
	start = clock();

	double beta, beta0, beta1, tolerance = 1e-4, gtheta, Fb0, Fb1, Fb, intm;
	MATRIX rbar0, rbar1, normal, rbar, pass;
	beta0 = 0.0;
	beta1 = 1.0;
	rbar0 = alpha_ns + beta0 * (r - alpha_ns);
	rbar1 = alpha_ns + beta1 * (r - alpha_ns);
	if ((r % r) < tolerance && (alpha_ns % alpha_ns) < tolerance)
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				normal(i, j) = 0.0;
		normal(0, 0) = -1 / sqrt(6.0);
		normal(1, 1) = -1.0 / sqrt(6.0);
		normal(2, 2) = 2.0 / sqrt(6.0);
		rbar = sqrt(2.0 / 3) * Mfc*exp(-np * psi)*normal;
		beta = 1.0e20;
	}
	else if (sqrt((r - alpha_ns) % (r - alpha_ns)) < tolerance)
	{
		normal = r1;
		rbar = sqrt(2.0 / 3) * etamplus1*sin3theta*normal;
		beta = 1.0e20;
	}
	else
	{
		if ((rbar0 % rbar0) < tolerance)
		{
			beta0 = 0.01;
			rbar0 = alpha_ns + beta0 * (r - alpha_ns);
		}
		normal = rbar0 / sqrt((rbar0 % rbar0));
		pass = normal * normal*normal;
		sin3theta = -sqrt(6.0)*tr(pass);
		if (sin3theta > 1.0)
			sin3theta = 1.0;
		else if (sin3theta < -1.0)
			sin3theta = -1.0;
		gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));
		Fb0 = (sqrt(2.0 / 3)*etamplus1*gtheta*normal - rbar0) % normal;
		normal = rbar1 / sqrt(rbar1 % rbar1);
		pass = normal * normal*normal;
		sin3theta = -sqrt(6.0)*tr(pass);
		if (sin3theta > 1.0)
			sin3theta = 1.0;
		else if (sin3theta < -1.0)
			sin3theta = -1.0;
		gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));
		Fb1 = (sqrt(2.0 / 3)*etamplus1*gtheta*normal - rbar1) % normal;
		if (abs(Fb0) <= 1.0e-5)
		{
			rbar = rbar0;
			beta = beta0;
		}
		else if (abs(Fb1) <= 1.0e-5)
		{
			rbar = rbar1;
			beta = beta1;
		}
		else
		{
			while (Fb0*Fb1 > 0)
			{
				beta0 = beta1;
				beta1 = 2 * beta1;
				rbar0 = alpha_ns + beta0 * (r - alpha_ns);
				rbar1 = alpha_ns + beta1 * (r - alpha_ns);
				normal = rbar0 / sqrt((rbar0 % rbar0));
				pass = normal * normal*normal;
				sin3theta = -sqrt(6.0)*tr(pass);
				if (sin3theta > 1.0)
					sin3theta = 1.0;
				else if (sin3theta < -1.0)
					sin3theta = -1.0;
				gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));
				Fb0 = (sqrt(2.0 / 3)*etamplus1*gtheta*normal - rbar0) % normal;
				normal = rbar1 / sqrt(rbar1 % rbar1);
				pass = normal * normal*normal;
				sin3theta = -sqrt(6.0)*tr(pass);
				if (sin3theta > 1.0)
					sin3theta = 1.0;
				else if (sin3theta < -1.0)
					sin3theta = -1.0;
				gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));
				Fb1 = (sqrt(2.0 / 3)*etamplus1*gtheta*normal - rbar1) % normal;
			}
			if (abs(Fb0) <= 1.0e-5)
			{
				rbar = rbar0;
				beta = beta0;
			}
			else if (abs(Fb1) <= 1.0e-5)
			{
				rbar = rbar1;
				beta = beta1;
			}
			else
			{
				beta = beta1 - Fb1 * (beta1 - beta0) / (Fb1 - Fb0);
				rbar = alpha_ns + beta * (r - alpha_ns);
				normal = rbar / sqrt(rbar % rbar);
				pass = normal * normal*normal;
				sin3theta = -sqrt(6.0)*tr(pass);
				if (sin3theta > 1.0)
					sin3theta = 1.0;
				else if (sin3theta < -1.0)
					sin3theta = -1.0;
				gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));
				Fb = (sqrt(2.0 / 3)*etamplus1*gtheta*normal - rbar) % normal;
				intm = 1;
				while (abs(Fb) > 1.0e-6)
				{
					if (Fb*Fb1 < 0)
					{
						beta0 = beta1;
						Fb0 = Fb1;
						beta1 = beta;
						Fb1 = Fb;
					}
					else
					{
						Fb0 = Fb1 * Fb0 / (Fb1 + Fb);
						beta1 = beta;
						Fb1 = Fb;
					}
					beta = beta1 - Fb1 * (beta1 - beta0) / (Fb1 - Fb0);
					rbar = alpha_ns + beta * (r - alpha_ns);
					normal = rbar / sqrt(rbar % rbar);
					pass = normal * normal*normal;
					sin3theta = -sqrt(6.0)*tr(pass);
					if (sin3theta > 1.0)
						sin3theta = 1.0;
					else if (sin3theta < -1.0)
						sin3theta = -1.0;
					gtheta = 1 / (1 + Mfc / 6.0*(sin3theta + sin3theta * sin3theta) + (Mfc - Mfo) / Mfo * (1 - sin3theta * sin3theta));
					Fb = (sqrt(2.0 / 3)*etamplus1*gtheta*normal - rbar) % normal;
					intm = intm + 1;
				}
			}
		}

	}
	
	finish = clock();
	betaTimer += (double)(finish - start) / CLOCKS_PER_SEC;
	return beta;
}
