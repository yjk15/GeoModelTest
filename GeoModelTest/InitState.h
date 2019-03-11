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

	QLabel *labelStress, *labelStrain;
	QLineEdit inputStrain[9], inputStress[9];
	QPushButton *yes, *cancel;
	MATRIX stress, strain;

signals:
	void sendInitState(MATRIX, MATRIX);

private slots:
	void clickYes();
	void clickCancel();
};