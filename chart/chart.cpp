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
	QValueAxis *axisX = new QValueAxis; //定义X轴
	axisX->setLabelFormat("%g"); //设置刻度的格式
	axisX->setTitleText("X Axis"); //设置X轴的标题
	axisX->setGridLineVisible(true); //设置是否显示网格线
	axisX->setMinorTickCount(4); //设置小刻度线的数目
   // axisX->setLabelsVisible(false); //设置刻度是否显示

	QValueAxis *axisY = new QValueAxis;
	axisY->setTitleText("Y Axis");
	axisY->setLabelFormat("%.2f");
	axisY->setGridLineVisible(true);

	chart->setAxisX(axisX, series);
	chart->setAxisY(axisY, series);
	chart->setTitle("nihao");

	widget->setChart(chart);
}
