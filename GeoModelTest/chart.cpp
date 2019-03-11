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
	for (int i = 0; i < 100; i++)
		series->append(i, sin(0.5 * i));
	chart = new QChart();
	chart->addSeries(series);
	chart->legend()->hide();
	axisX = new QValueAxis; //定义X轴
	axisX->setLabelFormat("%g"); //设置刻度的格式
	axisX->setTitleText("X Axis"); //设置X轴的标题
	axisX->setGridLineVisible(true); //设置是否显示网格线
	axisX->setMinorTickCount(4); //设置小刻度线的数目
   // axisX->setLabelsVisible(false); //设置刻度是否显示

	axisY = new QValueAxis;
	axisY->setTitleText("Y Axis");
	axisY->setLabelFormat("%.2f");
	axisY->setGridLineVisible(true);

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

	widget = new QChartView(this);
	widget->setGeometry(QRect(0, 0, 600, 400));
	chart = new QChart();
	series = new QSplineSeries();
	axisX = new QValueAxis; //定义X轴
	axisY = new QValueAxis;

	for (int i = 0; i < model->stressPath->size(); i++)
		series->append(getAxisX(model->axisX, i), getAxisY(model->axisY, i));

	chart->addSeries(series);
	chart->legend()->hide();
	axisX->setLabelFormat("%g"); //设置刻度的格式
	axisX->setTitleText(getAxisXTitle(model->axisX)); //设置X轴的标题
	axisX->setGridLineVisible(true); //设置是否显示网格线
	axisX->setMinorTickCount(4); //设置小刻度线的数目
	// axisX->setLabelsVisible(false); //设置刻度是否显示

	axisY->setTitleText(getAxisYTitle(model->axisY));
	axisY->setLabelFormat("%g");
	axisY->setGridLineVisible(true);
	axisY->setMinorTickCount(4);

	chart->setAxisX(axisX, series);
	chart->setAxisY(axisY, series);
	chart->setTitle(*(model->figureTitle));

	widget->setChart(chart);
}

double Chart::getAxisX(int axisX, int i) {
	switch (axisX) {
	case 0:
		return model->stressPath->at(i)(0, 0);
	case 1:
		return model->stressPath->at(i)(1, 1);
	case 2:
		return model->stressPath->at(i)(2, 2);
	case 3:
		return model->stressPath->at(i)(0, 0) - model->stressPath->at(i)(2, 2);
	default:
		return 0;
	}
}

double Chart::getAxisY(int axisY, int i) {
	switch (axisY) {
	case 0:
		return model->strainPath->at(i)(0, 0);
	case 1:
		return model->strainPath->at(i)(1, 1);
	case 2:
		return model->strainPath->at(i)(2, 2);
	case 3:
		return (model->strainPath->at(i)(0, 0) - model->strainPath->at(i)(2, 2)) * 2 / 3;
	default:
		return 0;
	}
}

QString Chart::getAxisXTitle(int axisX) {
	switch (axisX) {
	case 0:
		return "stress(0, 0)";
	case 1:
		return "stress(1, 1)";
	case 2:
		return "stress(2, 2)";
	case 3:
		return "q";
	default:
		return "";
	}
}

QString Chart::getAxisYTitle(int axisY) {
	switch (axisY) {
	case 0:
		return "strain(0, 0)";
	case 1:
		return "strain(1, 1)";
	case 2:
		return "strain(2, 2)";
	case 3:
		return "epsilonq";
	default:
		return "";
	}
}