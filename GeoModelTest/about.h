#pragma once

#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>

class ABOUT : public QDialog {
	Q_OBJECT
public:
	ABOUT(QWidget *parent = 0);
	~ABOUT();

	QLabel *label1;
	QLabel *label2;
	QLabel *label3;
};