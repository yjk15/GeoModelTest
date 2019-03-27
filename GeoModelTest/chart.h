#pragma once

#include <QtCharts>
#include "Model.h"
using namespace QtCharts;

#include <QtWidgets/QMainWindow>
#include "ui_chart.h"

class Chart : public QMainWindow
{
	Q_OBJECT

public:
	Chart(MODEL *m, QWidget *parent = Q_NULLPTR);
	~Chart();
	QChartView *widget;
	MODEL *model;
	QSplineSeries *series;
	QValueAxis *axisX;
	QChart *chart;
	QValueAxis *axisY;

public:
	void setNewChart();

private:
	Ui::chartClass ui;

private:
	double getAxis(int axis, int i);
	QString getAxisXTitle(int axisX);
	QString getAxisYTitle(int axisY);
};
