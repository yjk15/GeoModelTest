#pragma once

#include <QtCharts>
using namespace QtCharts;

#include <QtWidgets/QMainWindow>
#include "ui_chart.h"

class chart : public QMainWindow
{
	Q_OBJECT

public:
	chart(QWidget *parent = Q_NULLPTR);
	QChartView *widget;

private:
	Ui::chartClass ui;
};
