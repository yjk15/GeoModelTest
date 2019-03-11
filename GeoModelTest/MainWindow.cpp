#include "MainWindow.h"
#include <QtGui>
#include <QtWidgets>

#pragma execution_character_set("utf-8")    // ��������������⣬ע�⣡����

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	about = new ABOUT(this);
	inputModelE = new ModelParaInputE(this);
	initState = new InitState(this);
	endAndReversalState = new EndAndReversalState(this);
	endAndReversalStateNoCycle = new EndAndReversalStateNoCycle(this);
	model = new MODEL();
	figureTitle = new QLineEdit(this);
	figure = new Chart(model, this);
	displayParameter = new DisplayParameter(model, this);

	this->resize(500, 600);
	this->setWindowTitle("ɰ������ģ������");

	setMenu();
	ui.mainToolBar->hide();
	setType();
	setModel();
	calculate();
	drawFigure();
}

MainWindow::~MainWindow() {
	delete about, model, figureTitle;
}

void MainWindow::setMenu() {
	QMenu *menu = new QMenu("����");
	QAction *userManual = new QAction(this);
	userManual->setText("�û��ֲ�");
	menu->addAction(userManual);
	connect(userManual, SIGNAL(triggered()), this, SLOT(openManual()));

	QAction *about = new QAction(this);
	about->setText("����");
	menu->addAction(about);
	connect(about, SIGNAL(triggered()), this, SLOT(openAbout()));
	ui.menuBar->addMenu(menu);
}

void MainWindow::setType() {
	QPushButton *initStateButton = new QPushButton(this);
	initStateButton->setFont(QFont("Timers", 11, QFont::Courier));
	initStateButton->resize(80, 30);
	initStateButton->setText(tr("��ʼ״̬"));
	initStateButton->move(280, 70);
	connect(initStateButton, SIGNAL(clicked()), this, SLOT(setInitState()));

	QPushButton *endStateButton = new QPushButton(this);
	endStateButton->setFont(QFont("Timers", 11, QFont::Courier));
	endStateButton->resize(80, 30);
	endStateButton->setText(tr("ĩ״̬"));
	endStateButton->move(390, 70);
	connect(endStateButton, SIGNAL(clicked()), this, SLOT(setEndState()));

	QComboBox *testType = new QComboBox(this);
	testType->setFont(QFont("Timers", 11, QFont::Courier));
	testType->resize(180,30);
	testType->move(30, 70);
	testType->addItem(QWidget::tr("����ˮ����ѹ������"));
	testType->addItem(QWidget::tr("����ˮ���ἷ������"));
	testType->addItem(QWidget::tr("����ˮ����ѭ������"));
	testType->addItem(QWidget::tr("��ˮ����ѹ������"));
	testType->addItem(QWidget::tr("��ˮ���ἷ������"));
	testType->addItem(QWidget::tr("��ˮ����ѭ������"));
	connect(testType, SIGNAL(currentIndexChanged(int)), this, SLOT(setTestType(int)));
}

void MainWindow::setModel() {
	QComboBox *consModel = new QComboBox(this);
	consModel->setFont(QFont("Timers", 11, QFont::Courier));
	consModel->resize(180, 30);
	consModel->move(30, 164);
	consModel->height();
	consModel->addItem(QWidget::tr("����ģ��"));
	consModel->addItem(QWidget::tr("EBģ��"));
	consModel->addItem(QWidget::tr("Dafalias and Manzariģ��"));
	connect(consModel, SIGNAL(currentIndexChanged(int)), this, SLOT(setConsModel(int)));

	QPushButton *modelParameter = new QPushButton(this);
	modelParameter->setFont(QFont("Timers", 11, QFont::Courier));
	modelParameter->resize(80, 30);
	modelParameter->setText(tr("ģ�Ͳ���"));
	modelParameter->move(335, 164);
	connect(modelParameter, SIGNAL(clicked()), this, SLOT(setModelParameter()));
}

void MainWindow::calculate() {
	QLabel *labelStep= new QLabel(this);
	labelStep->setFont(QFont("Timers", 11, QFont::Courier));
	labelStep->setText("����");
	labelStep->resize(80, 30);
	labelStep->move(30, 265);

	QComboBox *step = new QComboBox(this);
	step->setFont(QFont("Timers", 11, QFont::Courier));
	step->resize(120, 30);
	step->move(80, 265);
	step->addItem("1e-5");
	step->addItem("5e-6");
	step->addItem("1e-6");
	step->addItem("5e-7");
	step->addItem("1e-7");
	connect(step, SIGNAL(currentIndexChanged(int)), this, SLOT(setStepLength(int)));

	QPushButton *check = new QPushButton(this);
	check->setFont(QFont("Timers", 11, QFont::Courier));
	check->resize(100, 30);
	check->setText(tr("��������"));
	check->move(80, 330);
	connect(check, SIGNAL(clicked()), this, SLOT(checkParameter()));

	QPushButton *calculateButton = new QPushButton(this);
	calculateButton->setFont(QFont("Timers", 11, QFont::Courier));
	calculateButton->resize(80, 80);
	calculateButton->setText(tr("��ʼ����"));
	calculateButton->move(335, 270);
	connect(calculateButton, SIGNAL(clicked()), this, SLOT(beginCalculate()));
}

void MainWindow::drawFigure() {
	figureTitle->resize(320, 30);
	figureTitle->move(110, 430);
	figureTitle->setMaxLength(30);
	figureTitle->setFont(QFont("Timers", 11, QFont::Courier));
	connect(figureTitle, SIGNAL(editingFinished()), this, SLOT(setFigureTitle()));

	QLabel *label1 = new QLabel(this);
	label1->setFont(QFont("Timers", 11, QFont::Courier));
	label1->setText("ͼ�����");
	label1->resize(80, 30);
	label1->move(30, 430);

	QLabel *label2 = new QLabel(this);
	label2->setFont(QFont("Timers", 11, QFont::Courier));
	label2->setText("X��");
	label2->resize(80, 30);
	label2->move(30, 480);

	QComboBox *axisX = new QComboBox(this);
	axisX->setFont(QFont("Timers", 11, QFont::Courier));
	axisX->resize(120, 30);
	axisX->move(80, 480);
	axisX->addItem("stress11");
	axisX->addItem("stress22");
	axisX->addItem("stress33");
	axisX->addItem("q");
	connect(axisX, SIGNAL(currentIndexChanged(int)), this, SLOT(setAxisX(int)));

	QLabel *label3 = new QLabel(this);
	label3->setFont(QFont("Timers", 11, QFont::Courier));
	label3->setText("Y��");
	label3->resize(80, 30);
	label3->move(285, 480);

	QComboBox *axisY = new QComboBox(this);
	axisY->setFont(QFont("Timers", 11, QFont::Courier));
	axisY->resize(120, 30);
	axisY->move(335, 480);
	axisY->addItem("strain11");
	axisY->addItem("strain22");
	axisY->addItem("strain33");
	axisY->addItem("epsilonq"); 
	connect(axisY, SIGNAL(currentIndexChanged(int)), this, SLOT(setAxisY(int)));

	QPushButton *draw = new QPushButton(this);
	draw->setFont(QFont("Timers", 11, QFont::Courier));
	draw->resize(80, 30);
	draw->setText(tr("��ͼ"));
	draw->move(80, 530);
	connect(draw, SIGNAL(clicked()), this, SLOT(drawing()));

	QPushButton *resetButton = new QPushButton(this);
	resetButton->setFont(QFont("Timers", 11, QFont::Courier));
	resetButton->resize(80, 30);
	resetButton->setText(tr("�������"));
	resetButton->move(335, 530);
	connect(resetButton, SIGNAL(clicked()), this, SLOT(reset()));
}

void MainWindow::setInitState() {
	ui.statusBar->showMessage(tr("������ʼ״̬"));
	connect(initState, SIGNAL(sendInitState(MATRIX, MATRIX)), this, SLOT(receiveInitState(MATRIX, MATRIX)));
	initState->show();
}

void MainWindow::receiveInitState(MATRIX stress, MATRIX strain) {
	model->stress = stress;
	model->strain = strain;
}

void MainWindow::setEndState() {
	ui.statusBar->showMessage(tr("������ֹ/��ת״̬"));
	if (model->testType == 2 || model->testType == 5) {
		connect(endAndReversalState, SIGNAL(sendEndAndReversalState(int, double, int)), this, SLOT(receiveEndAndReversalState(int, double, int)));
		endAndReversalState->show();
	}
	else {
		connect(endAndReversalStateNoCycle, SIGNAL(sendEndAndReversalState(int, double, int)), this, SLOT(receiveEndAndReversalState(int, double, int)));
		endAndReversalStateNoCycle->show();
	}
}

void MainWindow::receiveEndAndReversalState(int type, double point, int loop) {
	model->endAndReversalType = type;
	model->endAndReversalPoint = point;
	model->loop = loop;
}

void MainWindow::setTestType(int type) {
	ui.statusBar->showMessage(tr("������������"));
	this->model->testType = type;
	//char a[10];
	//_itoa(this->model->testType, a, 10);
	//QMessageBox::information(this, tr("hello"), tr(a));
}

void MainWindow::openManual() {
	QString qtManulFile = "./trans.pdf";
	QDesktopServices::openUrl(QUrl::fromLocalFile(qtManulFile));
}

void MainWindow::openAbout() {
	about->show();
}

void MainWindow::setConsModel(int model) {
	ui.statusBar->showMessage(tr("ѡ��ɰ������ģ��"));
	this->model->model = model;
}

void MainWindow::setModelParameter() {
	ui.statusBar->showMessage(tr("����ģ�Ͳ���"));
	switch (model->model) {
	case 0:
		connect(inputModelE, SIGNAL(sendEv(double, double)), this, SLOT(receiveE(double, double)));
		inputModelE->show();
		break;
	default:
		QMessageBox::information(this, tr("hello"), tr("Hello World!"));
	}
}

void MainWindow::receiveE(double E, double v) {
	model->internalParameter[0] = E;
	model->internalParameter[1] = v;
}

void MainWindow::setStepLength(int l) {
	ui.statusBar->showMessage(tr("���ò���"));
	switch (l) {
	case 0:
		this->model->stepLength = 1e-5;
		break;
	case 1:
		this->model->stepLength = 5e-6;
		break;
	case 2:
		this->model->stepLength = 1e-6;
		break;
	case 3:
		this->model->stepLength = 5e-7;
		break;
	case 4:
		this->model->stepLength = 1e-7;
		break;
	default:
		this->model->stepLength = 1e-7;
	}
	
}

void MainWindow::checkParameter() {
	ui.statusBar->showMessage(tr("������"));
	delete displayParameter;
	displayParameter = new DisplayParameter(model, this);
	displayParameter->show();
}

void MainWindow::beginCalculate() {
	ui.statusBar->showMessage(tr("�����У����������������"));
	model->Simulate();
	ui.statusBar->showMessage(tr("�������"));
}

void MainWindow::setFigureTitle() {
	*this->model->figureTitle = this->figureTitle->text();
}

void MainWindow::setAxisX(int x) {
	ui.statusBar->showMessage(tr("����X��"));
	this->model->axisX = x;
}

void MainWindow::setAxisY(int y) {
	ui.statusBar->showMessage(tr("����Y��"));
	this->model->axisY = y;
}

void MainWindow::drawing() {
	ui.statusBar->showMessage(tr("��ͼ"));
	figure->setNewChart();
	figure->show();
}

void MainWindow::reset() {
	delete model, figure, displayParameter;
	model = new MODEL();
	figure = new Chart(model, this);
	displayParameter = new DisplayParameter(model, this);
	ui.statusBar->showMessage(tr("���������в���"));
}