#pragma once
#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include "ui_ModelParaInputDMF.h"

class ModelParaInputDMF : public QDialog
{
	Q_OBJECT
public:
	ModelParaInputDMF(QWidget *parent = 0);
	~ModelParaInputDMF();


	QLineEdit *input;
	QPushButton *yes, *cancel;
	double para[17]; //“˛ Ω/œ‘ Ω, G0, v, M, c, lambdaC, e0, xi, m, h0, ch, nb, nd, A0, nd, Fin, cF, eA

private:
	Ui::DMF ui;

signals:
	void sendParaDMF(double *para);

private slots:
	void clickYes();
	void closeDialog();
};