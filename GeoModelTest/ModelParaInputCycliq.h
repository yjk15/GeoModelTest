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
	double para[15]; //解法（暂时只有一种）, G0, kappa, h, M, dre1, dre2, dir, alpha/eta, gammadr, np, nd, lambdaC, e0, xi
	//注：在论文中的alpha即为此程序中的eta
private:
	Ui::DialogCycliq ui;

signals:
	void sendParaCycliq(double *para);

private slots:
	void clickYes();
	void closeDialog();
};