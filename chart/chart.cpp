 #include "chart.h"

chart::chart(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	ui.centralWidget->resize(600, 400);
	ui.statusBar->hide();
	ui.mainToolBar->hide();
	ui.menuBar->hide();

	widget = new QChartView(this);
	widget->setGeometry(QRect(0, 0, 600, 400));

	QSplineSeries *series = new QSplineSeries();
	for (int i = 0; i < 100; i++)
		series->append(i, sin(0.5 * i));
	QChart *chart = new QChart();
	chart->addSeries(series);
	chart->legend()->hide();
	QValueAxis *axisX = new QValueAxis; //����X��
	axisX->setLabelFormat("%g"); //���ÿ̶ȵĸ�ʽ
	axisX->setTitleText("X Axis"); //����X��ı���
	axisX->setGridLineVisible(true); //�����Ƿ���ʾ������
	axisX->setMinorTickCount(4); //����С�̶��ߵ���Ŀ
   // axisX->setLabelsVisible(false); //���ÿ̶��Ƿ���ʾ

	QValueAxis *axisY = new QValueAxis;
	axisY->setTitleText("Y Axis");
	axisY->setLabelFormat("%.2f");
	axisY->setGridLineVisible(true);

	chart->setAxisX(axisX, series);
	chart->setAxisY(axisY, series);
	chart->setTitle("nihao");

	widget->setChart(chart);
}
