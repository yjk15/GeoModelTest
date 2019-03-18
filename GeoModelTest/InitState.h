#pragma once
#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
#include "matrix.h"

class InitState : public QDialog {
	Q_OBJECT
public:
	InitState(QWidget *parent = 0);
	~InitState();

	QLabel *labelStress, *labelStrain, *labelE;
	QLineEdit inputStrain[9], inputStress[9], *inputE;
	QPushButton *yes, *cancel;
	MATRIX stress, strain;
	double ee;

signals:
	void sendInitState(MATRIX, MATRIX, double);

private slots:
	void clickYes();
	void clickCancel();
};