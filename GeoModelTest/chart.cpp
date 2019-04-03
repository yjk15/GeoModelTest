 #include "chart.h"

#pragma execution_character_set("utf-8")    // ��������������⣬ע�⣡����

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
	axisX = new QValueAxis; //����X��
	axisX->setLabelFormat("%g"); //���ÿ̶ȵĸ�ʽ
	axisX->setTitleText("X Axis"); //����X��ı���
	axisX->setGridLineVisible(true); //�����Ƿ���ʾ������
	axisX->setMinorTickCount(4); //����С�̶��ߵ���Ŀ
   // axisX->setLabelsVisible(false); //���ÿ̶��Ƿ���ʾ

	axisY = new QValueAxis;
	axisY->setTitleText("Y Axis");
	axisY->setLabelFormat("%g");
	axisY->setGridLineVisible(true);
	axisY->setMinorTickCount(4); //����С�̶��ߵ���Ŀ

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
	series = new QSplineSeries(this);
	axisX = new QValueAxis(this); //����X��
	axisY = new QValueAxis(this);

	for (int i = 0; i < model->stressPath->size(); i++)
		series->append(getAxis(model->axisX, i), getAxis(model->axisY, i));

	chart->addSeries(series);
	chart->legend()->hide();
	axisX->setLabelFormat("%g"); //���ÿ̶ȵĸ�ʽ
	axisX->setTitleText(getAxisTitle(model->axisX)); //����X��ı���
	axisX->setGridLineVisible(true); //�����Ƿ���ʾ������
	axisX->setMinorTickCount(4); //����С�̶��ߵ���Ŀ
	// axisX->setLabelsVisible(false); //���ÿ̶��Ƿ���ʾ

	axisY->setTitleText(getAxisTitle(model->axisY));
	axisY->setLabelFormat("%g");
	axisY->setGridLineVisible(true);
	axisY->setMinorTickCount(4);

	chart->setAxisX(axisX, series);
	chart->setAxisY(axisY, series);
	chart->setTitle(*(model->figureTitle));

	widget->setChart(chart);
}

double Chart::getAxis(int axis, int i) {
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
		return "��q";
	case 11:
		return "��p";
	case 12:
		return "e";
	default:
		return "";
	}
}