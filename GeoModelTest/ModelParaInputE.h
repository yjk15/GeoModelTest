#pragma once
#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>

class ModelParaInputE : public QDialog
{
	Q_OBJECT
public:
	ModelParaInputE(QWidget *parent = 0);
	~ModelParaInputE();

	QLabel *labelE, *labelv, *labelGPa;
	QLineEdit *inputE, *inputv;
	QPushButton *yes, *cancel;
	QString E, v;

signals:
	void sendEv(double E, double v);

private slots:
	void clickYes();
	void closeDialog();
};