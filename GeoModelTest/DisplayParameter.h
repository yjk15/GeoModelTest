#pragma once

#include <QWidget>
#include <QDebug>
#include <QPushButton>
#include <QDialog>
#include <QLabel>
#include <QString>
#include "Model.h"

class DisplayParameter : public QDialog {
	Q_OBJECT
public:
	DisplayParameter(MODEL *m, QWidget *parent = 0);
	~DisplayParameter();

	MODEL *model;

	QLabel *label, labelee[4];

	void DisplayPara(int modelType, int testType);
};