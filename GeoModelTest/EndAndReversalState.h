#pragma once
#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>

class EndAndReversalState : public QDialog {
	Q_OBJECT
public:
	EndAndReversalState(QWidget *parent = 0);
	~EndAndReversalState();

	QLabel *labelType, *labelPoint, *labelLoop;
	QComboBox *inputType;
	QLineEdit *inputPoint, *inputLoop;
	QPushButton *yes, *cancel;

signals:
	void sendEndAndReversalState(int, double, int);

private slots:
	void clickYes();
	void clickCancel();
};