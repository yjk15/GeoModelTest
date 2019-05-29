 #include "chart.h"

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

Chart::Chart(MODEL *m, QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	model = m;

	ui.centralWidget->resize(600, 400);
	ui.statusBar->hide();
	ui.mainToolBar->hide();
	ui.menuBar->hide();

	widget = new QChartView(this);
	widget->setGeometry(QRect(0, 0, 600, 400));

	series = new QSplineSeries();
	for (int i = 0; i < model->stressPath->size(); i++)
		series->append(getAxis(model->axisX, i), getAxis(model->axisY, i));
	chart = new QChart();
	chart->addSeries(series);
	chart->legend()->hide();
	axisX = new QValueAxis; //定义X轴
	axisX->setLabelFormat("%g"); //设置刻度的格式
	axisX->setTitleText(getAxisTitle(model->axisX)); //设置X轴的标题
	axisX->setGridLineVisible(false); //设置是否显示网格线
	//axisX->setMinorTickCount(4); //设置小刻度线的数目
   // axisX->setLabelsVisible(false); //设置刻度是否显示

	axisY = new QValueAxis;
	axisY->setTitleText(getAxisTitle(model->axisY));
	axisY->setLabelFormat("%g");
	axisY->setGridLineVisible(false);
	//axisY->setMinorTickCount(4); //设置小刻度线的数目

	chart->setAxisX(axisX, series);
	chart->setAxisY(axisY, series);
	chart->setTitle(*(model->figureTitle));

	widget->setChart(chart);
}

Chart::~Chart() {
	delete series, chart, axisX, axisY, widget;
}

void Chart::setNewChart(){
	delete series, chart, axisX, axisY, widget;

	double minX, maxX, minY, maxY;

	widget = new QChartView(this);
	widget->setGeometry(QRect(0, 0, 600, 400));
	chart = new QChart();
	series = new QSplineSeries(this);
	axisX = new QValueAxis(this); //定义X轴
	axisY = new QValueAxis(this);

	minX = getAxis(model->axisX, 0);
	maxX = getAxis(model->axisX, 0);
	minY = getAxis(model->axisY, 0);
	maxY = getAxis(model->axisY, 0);
	for (int i = 0; i < model->stressPath->size(); i++) {
		series->append(getAxis(model->axisX, i), getAxis(model->axisY, i));
		//if (getAxis(model->axisX, i) < minX)
		//	minX = getAxis(model->axisX, i);
		//if (getAxis(model->axisX, i) > maxX)
		//	maxX = getAxis(model->axisX, i);
		//if (getAxis(model->axisY, i) < minY)
		//	minY = getAxis(model->axisY, i);
		//if (getAxis(model->axisY, i) > maxY)
		//	maxY = getAxis(model->axisY, i);
	}

	/*int E1 = 0;
	while (abs(minX) < 100) {
		minX *= 10;
		E1 += 1;
	}
	minX = minX / pow(10, E1);
	double E2 = 0;
	while (abs(maxX) < 100) {
		maxX *= 10;
		E2 += 1;
	}
	maxX = maxX / pow(10, E2);
	if (E2 < E1)
		E1 = E2;
	minX = double(int(minX * pow(10, E1) / 5) * 5) / pow(10, E1);
	maxX = double(int(maxX * pow(10, E1) / 5) * 5) / pow(10, E1);

	E1 = 0;
	while (abs(minY) < 100) {
		minY *= 10;
		E1 += 1;
	}
	minY = minY / pow(10, E1);
	E2 = 0;
	while (abs(maxY) < 100) {
		maxY *= 10;
		E2 += 1;
	}
	maxY = maxY / pow(10, E2);
	if (E2 < E1)
		E1 = E2;
	minY = double(int(minY * pow(10, E1) / 5) * 5) / pow(10, E1);
	maxY = double(int(maxY * pow(10, E1) / 5) * 5) / pow(10, E1);*/

	chart->addSeries(series);
	chart->legend()->hide();
	axisX->setLabelFormat("%g"); //设置刻度的格式
	axisX->setTitleText(getAxisTitle(model->axisX)); //设置X轴的标题
	axisX->setGridLineVisible(true); //设置是否显示网格线
	//axisX->setMinorTickCount(4); //设置小刻度线的数目
	// axisX->setLabelsVisible(false); //设置刻度是否显示
	/*if(model->axisX == 5)
		axisX->setRange(0, 100);*/
	//axisX->setRange(minX, maxX);

	axisY->setTitleText(getAxisTitle(model->axisY));
	axisY->setLabelFormat("%g");
	axisY->setGridLineVisible(true);
	//axisY->setMinorTickCount(4);
	/*if (model->axisY == 4)
		axisY->setRange(-40, 40);*/
	//axisY->setRange(minY, maxY);

	chart->setAxisX(axisX, series);
	chart->setAxisY(axisY, series);
	chart->setTitle(*(model->figureTitle));

	widget->setChart(chart);
}

double Chart::getAxis(int axis, int i) {
	MATRIX I(1,1,1);
	double p;
	switch (axis) {
	case 0:
		return model->stressPath->at(i)(0, 0);
	case 1:
		return model->stressPath->at(i)(1, 1);
	case 2:
		return model->stressPath->at(i)(2, 2);
	case 3:
		return model->stressPath->at(i)(0, 2);
	case 4:
		return model->stressPath->at(i)(2, 2) - model->stressPath->at(i)(0, 0);
	case 5:
		return (model->stressPath->at(i)(0, 0) + model->stressPath->at(i)(1, 1) + model->stressPath->at(i)(2, 2)) / 3.0;
	case 6:
		return model->strainPath->at(i)(0, 0);
	case 7:
		return model->strainPath->at(i)(1, 1);
	case 8:
		return model->strainPath->at(i)(2, 2);
	case 9:
		return model->strainPath->at(i)(0, 2);
	case 10:
		return (model->strainPath->at(i)(2, 2) - model->strainPath->at(i)(0, 0)) * 2.0 / 3;
	case 11:
		return model->strainPath->at(i)(0, 0) + model->strainPath->at(i)(1, 1) + model->strainPath->at(i)(2, 2);
	case 12:
		return model->saveParameter.at(i).at(0);
	default:
		return 0;
	}
}

QString Chart::getAxisTitle(int axisX) {
	switch (axisX) {
	case 0:
		return "stress(0, 0)";
	case 1:
		return "stress(1, 1)";
	case 2:
		return "stress(2, 2)";
	case 3:
		return "stress(1, 3)";
	case 4:
		return "q";
	case 5:
		return "p";
	case 6:
		return "strain(0, 0)";
	case 7:
		return "strain(1, 1)";
	case 8:
		return "strain(2, 2)";
	case 9:
		return "strain(1, 3)";
	case 10:
		return "εq";
	case 11:
		return "εp";
	case 12:
		return "e";
	default:
		return "";
	}
}