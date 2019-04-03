#pragma once
#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include "ui_ModelParaInputCycliq.h"

class ModelParaInputCycliq : public QDialog
{
	Q_OBJECT
public:
	ModelParaInputCycliq(QWidget *parent = 0);
	~ModelParaInputCycliq();


	QLineEdit *input;
	QComboBox *method;
	QPushButton *yes, *cancel;
	double para[15]; //�ⷨ����ʱֻ��һ�֣�, G0, kappa, h, M, dre1, dre2, dir, alpha/eta, gammadr, np, nd, lambdaC, e0, xi
	//ע���������е�alpha��Ϊ�˳����е�eta
private:
	Ui::DialogCycliq ui;

signals:
	void sendParaCycliq(double *para);

private slots:
	void clickYes();
	void closeDialog();
};