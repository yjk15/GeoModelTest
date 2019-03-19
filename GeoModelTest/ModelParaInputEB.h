#pragma once
#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>

class ModelParaInputEB : public QDialog
{
	Q_OBJECT
public:
	ModelParaInputEB(QWidget *parent = 0);
	~ModelParaInputEB();

	QLabel *label0, *label1, *label2, *label3, *label4;
	QLineEdit *inputEB0, *inputEB1, *inputEB2, *inputEB3, *inputEB4;
	QPushButton *yes, *cancel;
	QString Kd1, Kd2, nd2, nud, gammamax;

signals:
	void sendEBPara(double Kd1, double Kd2, double nd2, double nud, double gammamax);

private slots:
	void clickYes();
	void closeDialog();
};