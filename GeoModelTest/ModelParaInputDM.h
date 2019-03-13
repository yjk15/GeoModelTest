#pragma once
#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include "ui_ModelParaInputDM.h"

class ModelParaInputDM : public QDialog
{
	Q_OBJECT
public:
	ModelParaInputDM(QWidget *parent = 0);
	~ModelParaInputDM();

	
	QLineEdit *input;
	QComboBox *method;
	QPushButton *yes, *cancel;
	double para[16]; //“˛ Ω/œ‘ Ω, G0, v, M, c, lambdaC, e0, xi, m, h0, ch, nb, nd, A0, nd, zmax, cz

private:
	Ui::Dialog ui;

signals:
	void sendParaDM(double *para);

private slots:
	void clickYes();
	void closeDialog();
};