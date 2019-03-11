#pragma once
#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>

class EndAndReversalStateNoCycle : public QDialog {
	Q_OBJECT
public:
	EndAndReversalStateNoCycle(QWidget *parent = 0);
	~EndAndReversalStateNoCycle();

	QLabel *labelType, *labelPoint;
	QComboBox *inputType;
	QLineEdit *inputPoint;
	QPushButton *yes, *cancel;

signals:
	void sendEndAndReversalState(int, double, int);

private slots:
	void clickYes();
	void clickCancel();
};