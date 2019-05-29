#include "MainWindow.h"
#include <QtGui>
#include <QtWidgets>

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	about = new ABOUT(this);
	filePath = "";
	resultFilePath = "";
	inputModelE = new ModelParaInputE(this);
	inputModelDM = new ModelParaInputDM(this);
	inputModelDMF = new ModelParaInputDMF(this);
	inputModelEB = new ModelParaInputEB(this);
	inputModelCycliq = new ModelParaInputCycliq(this);
	initState = new InitState(this);
	endAndReversalState = new EndAndReversalState(this);
	endAndReversalStateNoCycle = new EndAndReversalStateNoCycle(this);
	model = new MODEL();
	figureTitle = new QLineEdit(this);
	figure = new Chart(model, this);
	displayParameter = new DisplayParameter(model, this);

	this->resize(500, 600);
	this->setWindowTitle("GeoModelTest");
	this->setWindowIcon(QIcon("./GMT.ico"));

	setMenu();
	ui.mainToolBar->hide();
	setType();
	setModel();
	calculate();
	drawFigure();
}

MainWindow::~MainWindow() {
	delete about, model, figureTitle, inputModelE, inputModelDM, inputModelEB;
}

void MainWindow::setMenu() {
	QMenu *file = new QMenu("file");
	QAction *loadPara = new QAction(this);
	loadPara->setText("load parameter");
	file->addAction(loadPara);
	connect(loadPara, SIGNAL(triggered()), this, SLOT(loadParameter()));

	QAction *savePara = new QAction(this);
	savePara->setText("save parameter");
	file->addAction(savePara);
	connect(savePara, SIGNAL(triggered()), this, SLOT(saveParameter()));

	QAction *saveParaNew = new QAction(this);
	saveParaNew->setText("save parameter...");
	file->addAction(saveParaNew);
	connect(saveParaNew, SIGNAL(triggered()), this, SLOT(saveParameterInNewFile()));

	QAction *saveOutput = new QAction(this);
	saveOutput->setText("save result");
	file->addAction(saveOutput);
	connect(saveOutput, SIGNAL(triggered()), this, SLOT(saveResult()));
	
	ui.menuBar->addMenu(file);

	QMenu *menu = new QMenu("help");
	QAction *userManual = new QAction(this);
	userManual->setText("user manual");
	menu->addAction(userManual);
	connect(userManual, SIGNAL(triggered()), this, SLOT(openManual()));

	QAction *about = new QAction(this);
	about->setText("about");
	menu->addAction(about);
	connect(about, SIGNAL(triggered()), this, SLOT(openAbout()));

	QAction *test = new QAction(this);
	test->setText("test");
	menu->addAction(test);
	connect(test, SIGNAL(triggered()), this, SLOT(test()));
	
	ui.menuBar->addMenu(menu);
}

void MainWindow::setType() {
	QPushButton *initStateButton = new QPushButton(this);
	initStateButton->setFont(QFont("Timers", 11, QFont::Courier));
	initStateButton->resize(80, 30);
	initStateButton->setText(tr("initial"));
	initStateButton->move(335, 50);
	connect(initStateButton, SIGNAL(clicked()), this, SLOT(setInitState()));

	QPushButton *endStateButton = new QPushButton(this);
	endStateButton->setFont(QFont("Timers", 11, QFont::Courier));
	endStateButton->resize(80, 30);
	endStateButton->setText(tr("end/rev"));
	endStateButton->move(335, 90);
	connect(endStateButton, SIGNAL(clicked()), this, SLOT(setEndState()));

	QComboBox *testType = new QComboBox(this);
	testType->setFont(QFont("Timers", 11, QFont::Courier));
	testType->resize(280,30);
	testType->move(30, 70);
	//testType->addItem(QWidget::tr("不排水三轴压缩试验"));
	//testType->addItem(QWidget::tr("不排水三轴挤长试验"));
	//testType->addItem(QWidget::tr("不排水三轴循环试验"));
	//testType->addItem(QWidget::tr("排水三轴压缩试验"));
	//testType->addItem(QWidget::tr("排水三轴挤长试验"));
	//testType->addItem(QWidget::tr("不排水循环扭剪试验"));
	testType->addItem(QWidget::tr("undrained triaxial compression"));
	testType->addItem(QWidget::tr("undrained triaxial extension"));
	testType->addItem(QWidget::tr("cyclic undrained triaxial"));
	testType->addItem(QWidget::tr("drained triaxial compression"));
	testType->addItem(QWidget::tr("drained triaxial extension"));
	testType->addItem(QWidget::tr("undrained cyclic torsional"));
	connect(testType, SIGNAL(currentIndexChanged(int)), this, SLOT(setTestType(int)));
}

void MainWindow::setModel() {
	QComboBox *consModel = new QComboBox(this);
	consModel->setFont(QFont("Timers", 11, QFont::Courier));
	consModel->resize(280, 30);
	consModel->move(30, 164);
	consModel->height();
	consModel->addItem(QWidget::tr("linear elastic"));
	consModel->addItem(QWidget::tr("EB"));
	consModel->addItem(QWidget::tr("Dafalias and Manzari"));
	consModel->addItem(QWidget::tr("Cycliq"));
	consModel->addItem(QWidget::tr("DM with fabric"));
	connect(consModel, SIGNAL(currentIndexChanged(int)), this, SLOT(setConsModel(int)));

	QPushButton *modelParameter = new QPushButton(this);
	modelParameter->setFont(QFont("Timers", 11, QFont::Courier));
	modelParameter->resize(80, 30);
	modelParameter->setText(tr("parameter"));
	modelParameter->move(335, 164);
	connect(modelParameter, SIGNAL(clicked()), this, SLOT(setModelParameter()));
}

void MainWindow::calculate() {
	QLabel *labelStep= new QLabel(this);
	labelStep->setFont(QFont("Timers", 11, QFont::Courier));
	labelStep->setText("step length");
	labelStep->resize(120, 30);
	labelStep->move(30, 265);

	QComboBox *step = new QComboBox(this);
	step->setFont(QFont("Timers", 11, QFont::Courier));
	step->resize(120, 30);
	step->move(140, 265);
	step->addItem("1e-3");
	step->addItem("1e-4");
	step->addItem("1e-5");
	step->addItem("5e-6");
	step->addItem("1e-6");
	step->addItem("5e-7");
	step->addItem("1e-7");
	connect(step, SIGNAL(currentIndexChanged(int)), this, SLOT(setStepLength(int)));

	QPushButton *check = new QPushButton(this);
	check->setFont(QFont("Timers", 11, QFont::Courier));
	check->resize(100, 30);
	check->setText(tr("check"));
	check->move(140, 330);
	connect(check, SIGNAL(clicked()), this, SLOT(checkParameter()));

	QPushButton *calculateButton = new QPushButton(this);
	calculateButton->setFont(QFont("Timers", 11, QFont::Courier));
	calculateButton->resize(80, 80);
	calculateButton->setText(tr("calculate"));
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
	label1->setText("title");
	label1->resize(80, 30);
	label1->move(30, 430);

	QLabel *label2 = new QLabel(this);
	label2->setFont(QFont("Timers", 11, QFont::Courier));
	label2->setText("X轴");
	label2->resize(80, 30);
	label2->move(30, 480);

	QComboBox *axisX = new QComboBox(this);
	axisX->setFont(QFont("Timers", 11, QFont::Courier));
	axisX->resize(120, 30);
	axisX->move(80, 480);
	axisX->addItem("stress11");
	axisX->addItem("stress22");
	axisX->addItem("stress33");
	axisX->addItem("stress13");
	axisX->addItem("q");
	axisX->addItem("p");
	axisX->addItem("strain11");
	axisX->addItem("strain22");
	axisX->addItem("strain33");
	axisX->addItem("strain13");
	axisX->addItem("εq");
	axisX->addItem("εv");
	axisX->addItem("e");
	connect(axisX, SIGNAL(currentIndexChanged(int)), this, SLOT(setAxisX(int)));

	QLabel *label3 = new QLabel(this);
	label3->setFont(QFont("Timers", 11, QFont::Courier));
	label3->setText("Y轴");
	label3->resize(80, 30);
	label3->move(285, 480);

	QComboBox *axisY = new QComboBox(this);
	axisY->setFont(QFont("Timers", 11, QFont::Courier));
	axisY->resize(120, 30);
	axisY->move(335, 480);
	axisY->addItem("stress11");
	axisY->addItem("stress22");
	axisY->addItem("stress33");
	axisY->addItem("stress13");
	axisY->addItem("q");
	axisY->addItem("p");
	axisY->addItem("strain11");
	axisY->addItem("strain22");
	axisY->addItem("strain33");
	axisY->addItem("strain13");
	axisY->addItem("εq"); 
	axisY->addItem("εv");
	axisY->addItem("e");
	connect(axisY, SIGNAL(currentIndexChanged(int)), this, SLOT(setAxisY(int)));

	QPushButton *draw = new QPushButton(this);
	draw->setFont(QFont("Timers", 11, QFont::Courier));
	draw->resize(80, 30);
	draw->setText(tr("figure"));
	draw->move(80, 530);
	connect(draw, SIGNAL(clicked()), this, SLOT(drawing()));

	QPushButton *resetButton = new QPushButton(this);
	resetButton->setFont(QFont("Timers", 11, QFont::Courier));
	resetButton->resize(80, 30);
	resetButton->setText(tr("reset"));
	resetButton->move(335, 530);
	connect(resetButton, SIGNAL(clicked()), this, SLOT(reset()));
}

void MainWindow::setInitState() {
	ui.statusBar->showMessage(tr("set initial state"));
	connect(initState, SIGNAL(sendInitState(MATRIX, MATRIX, double)), this, SLOT(receiveInitState(MATRIX, MATRIX, double)));
	initState->show();
}

void MainWindow::receiveInitState(MATRIX stress, MATRIX strain, double ee) {
	model->stress = stress;
	model->strain = strain;
	model->ee = ee;
}

void MainWindow::setEndState() {
	//ui.statusBar->showMessage(tr("设置终止/反转状态"));
	ui.statusBar->showMessage(tr("set reversal/ending state"));
	if (model->testType == 2 || model->testType == 5) {
		connect(endAndReversalState, SIGNAL(sendEndAndReversalState(int, double, int)), this, SLOT(receiveEndAndReversalState(int, double, int)));
		endAndReversalState->show();
	}
	else {
		connect(endAndReversalStateNoCycle, SIGNAL(sendEndAndReversalState(int, double, int)), this, SLOT(receiveEndAndReversalState(int, double, int)));
		endAndReversalStateNoCycle->show();
	}
}

void MainWindow::receiveEndAndReversalState(int type, double point, int reverse) {
	model->endAndReversalType = type;
	model->endAndReversalPoint = point;
	model->reverse = reverse;
}

void MainWindow::setTestType(int type) {
	ui.statusBar->showMessage(tr("set test type"));
	this->model->testType = type;
	//char a[10];
	//_itoa(this->model->testType, a, 10);
	//QMessageBox::information(this, tr("hello"), tr(a));
}

void MainWindow::openManual() {
	QString qtManulFile = "./manual.pdf";
	QDesktopServices::openUrl(QUrl::fromLocalFile(qtManulFile));
}

void MainWindow::openAbout() {
	about->show();
}

void MainWindow::loadParameter() {
	filePath = QFileDialog::getOpenFileName(this, tr("load"), "", tr("data(*.dat);;txt(*.txt);;all(*.*)"));
	loadPara();
	ui.statusBar->showMessage(tr("parameter has been loaded"));
}

void MainWindow::loadPara() {
	QFile file(filePath);
	QString para;
	file.open(QIODevice::ReadOnly);
	para = file.readLine();
	model->testType = para.toInt();
	para = file.readLine();
	model->model = para.toInt();
	para = file.readLine();
	model->ee = para.toDouble();
	para = file.readLine();
	model->endAndReversalType = para.toInt();
	para = file.readLine();
	model->endAndReversalPoint = para.toDouble();
	para = file.readLine();
	model->reverse = para.toInt();
	para = file.readLine();
	model->reverseCounter = para.toInt();
	para = file.readLine();
	model->stepCounter = para.toInt();
	para = file.readLine();
	QVariant tmp = para;
	model->direction = tmp.toBool();
	para = file.readLine();
	model->axisX = para.toInt();
	para = file.readLine();
	model->axisY = para.toInt();
	para = file.readLine();
	model->stepLength = para.toDouble();
	for (int i = 0; i < 30; i++) {
		para = file.readLine();
		model->internalParameter[i] = para.toDouble();
	}
	for (int i = 0; i < 9; i++) {
		para = file.readLine();
		model->stress.matrix[i] = para.toDouble();
	}
	for (int i = 0; i < 9; i++) {
		para = file.readLine();
		model->strain.matrix[i] = para.toDouble();
	}
	para = file.readLine();
	*model->figureTitle = para;
	file.close();
}

void MainWindow::saveParameter() {
	if (filePath == "")
		filePath = QFileDialog::getSaveFileName(this, tr("save"), "", tr("data(*.dat);;txt(*.txt);;all(*.*)"));
	QString para;
	para.append(QString::number(model->testType));
	para.append("\n");
	para.append(QString::number(model->model));
	para.append("\n");
	para.append(QString::number(model->ee));
	para.append("\n");
	para.append(QString::number(model->endAndReversalType));
	para.append("\n");
	para.append(QString::number(model->endAndReversalPoint));
	para.append("\n");
	para.append(QString::number(model->reverse));
	para.append("\n");
	para.append(QString::number(0));
	para.append("\n");
	para.append(QString::number(0));
	para.append("\n");
	para.append(QString::number(model->direction));
	para.append("\n");
	para.append(QString::number(model->axisX));
	para.append("\n");
	para.append(QString::number(model->axisY));
	para.append("\n");
	para.append(QString::number(model->stepLength));
	para.append("\n");
	for (int i = 0; i < 30; i++) {
		para.append(QString::number(model->internalParameter[i]));
		para.append("\n");
	}
	for (int i = 0; i < 9; i++) {
		para.append(QString::number(model->stress.matrix[i]));
		para.append("\n");
	}
	for (int i = 0; i < 9; i++) {
		para.append(QString::number(model->strain.matrix[i]));
		para.append("\n");
	}
	para.append(*model->figureTitle);

	QFile file(filePath);
	file.open(QIODevice::WriteOnly);
	file.write(para.toUtf8());
	file.close();
	ui.statusBar->showMessage(tr("paramenter has been saved"));
}

void MainWindow::saveParameterInNewFile() {
	filePath = QFileDialog::getSaveFileName(this, tr("save"), "", tr("data(*.dat);;txt(*.txt);;all(*.*)"));
	QString para;
	para.append(QString::number(model->testType));
	para.append("\n");
	para.append(QString::number(model->model));
	para.append("\n");
	para.append(QString::number(model->ee));
	para.append("\n");
	para.append(QString::number(model->endAndReversalType));
	para.append("\n");
	para.append(QString::number(model->endAndReversalPoint));
	para.append("\n");
	para.append(QString::number(model->reverse));
	para.append("\n");
	para.append(QString::number(0));
	para.append("\n");
	para.append(QString::number(0));
	para.append("\n");
	para.append(QString::number(model->direction));
	para.append("\n");
	para.append(QString::number(model->axisX));
	para.append("\n");
	para.append(QString::number(model->axisY));
	para.append("\n");
	para.append(QString::number(model->stepLength));
	para.append("\n");
	for (int i = 0; i < 30; i++) {
		para.append(QString::number(model->internalParameter[i]));
		para.append("\n");
	}
	for (int i = 0; i < 9; i++) {
		para.append(QString::number(model->stress.matrix[i]));
		para.append("\n");
	}
	for (int i = 0; i < 9; i++) {
		para.append(QString::number(model->strain.matrix[i]));
		para.append("\n");
	}
	para.append(*model->figureTitle);

	QFile file(filePath);
	file.open(QIODevice::WriteOnly);
	file.write(para.toUtf8());
	file.close();
	//ui.statusBar->showMessage(tr("参数已存储于新文件"));
	ui.statusBar->showMessage(tr("parameter has been saved in a new document"));
}

void MainWindow::saveResult() {
	resultFilePath = QFileDialog::getSaveFileName(this, tr("save"), "", tr("data(*.dat);;txt(*.txt);;all(*.*)"));
	saveRes();
	//ui.statusBar->showMessage(tr("计算结果已存储"));
	ui.statusBar->showMessage(tr("result is saved"));
}

void MainWindow::saveRes() {
	QFile file(resultFilePath);
	file.open(QIODevice::WriteOnly);
	QString result;
	for (int i = 0; i < model->stressPath->size(); i++) {
		for (int j = 0; j < 9; j++) {
			result.append(QString::number(model->stressPath->at(i).matrix[j]));
			result.append(" ");
		}
		for (int j = 0; j < 9; j++) {
			result.append(QString::number(model->strainPath->at(i).matrix[j]));
			result.append(" ");
		}
		for (int j = 0; j < 2; j++) {
			result.append(QString::number(model->saveParameter.at(i).at(j)));
			result.append(" ");
		}
		result.append("\n");
	}
	file.write(result.toUtf8());
	file.close();
}

void MainWindow::setConsModel(int model) {
	ui.statusBar->showMessage(tr("选择砂土本构模型"));
	this->model->model = model;
}

void MainWindow::setModelParameter() {
	//ui.statusBar->showMessage(tr("设置模型参数"));
	ui.statusBar->showMessage(tr("set model parameter"));
	switch (model->model) {
	case 0:
		connect(inputModelE, SIGNAL(sendEv(double, double)), this, SLOT(receiveE(double, double)));
		inputModelE->show();
		break;
	case 1:
		connect(inputModelEB, SIGNAL(sendEBPara(double, double, double, double, double)), this, SLOT(receiveEB(double, double, double, double, double)));
		inputModelEB->show();
		break;
	case 2:
		connect(inputModelDM, SIGNAL(sendParaDM(double*)), this, SLOT(receiveDM(double*)));
		inputModelDM->show();
		break;
	case 3:
		connect(inputModelCycliq, SIGNAL(sendParaCycliq(double*)), this, SLOT(receiveCycliq(double*)));
		inputModelCycliq->show();
		break;
	case 4:
		connect(inputModelCycliq, SIGNAL(sendParaDMF(double*)), this, SLOT(receiveDMF(double*)));
		inputModelDMF->show();
		break;
	default:
		QMessageBox::information(this, tr("hello"), tr("Hello World!"));
	}
}

void MainWindow::receiveE(double E, double v) {
	model->internalParameter[0] = E;
	model->internalParameter[1] = v;
}

void MainWindow::receiveEB(double Kd1, double Kd2, double nd2, double nud, double gammamax) {
	model->internalParameter[0] = Kd1;
	model->internalParameter[1] = Kd2;
	model->internalParameter[2] = nd2;
	model->internalParameter[3] = nud;
	model->internalParameter[4] = gammamax;
}

void MainWindow::receiveDM(double* para) {
	for (int i = 0; i < 16; i++) {
		model->internalParameter[i] = para[i];
	}
}

void MainWindow::receiveDMF(double* para) {
	for (int i = 0; i < 17; i++) {
		model->internalParameter[i] = para[i];
	}
}

void MainWindow::receiveCycliq(double* para) {
	for (int i = 0; i < 15; i++) {
		model->internalParameter[i] = para[i];
	}
}

void MainWindow::setStepLength(int l) {
	//ui.statusBar->showMessage(tr("设置步长"));
	ui.statusBar->showMessage(tr("set step length"));
	switch (l) {
	case 0:
		this->model->stepLength = 1e-3;
		break;
	case 1:
		this->model->stepLength = 1e-4;
		break;
	case 2:
		this->model->stepLength = 1e-5;
		break;
	case 3:
		this->model->stepLength = 5e-6;
		break;
	case 4:
		this->model->stepLength = 1e-6;
		break;
	case 5:
		this->model->stepLength = 5e-7;
		break;
	case 6:
		this->model->stepLength = 1e-7;
		break;
	default:
		this->model->stepLength = 1e-7;
	}
	
}

void MainWindow::checkParameter() {
	//ui.statusBar->showMessage(tr("检查参数"));
	ui.statusBar->showMessage(tr("check parameter"));
	delete displayParameter;
	displayParameter = new DisplayParameter(model, this);
	displayParameter->show();
}

void MainWindow::beginCalculate() {
	QMessageBox *message = new QMessageBox(QMessageBox::Information, "calculating", tr("calculating..."), QMessageBox::NoButton, this);
	message->show();
	clock_t start, finish;
	double cost;
	QString tmp, out;
	start = clock();
	model->Simulate();
	finish = clock();
	cost = (double)(finish - start) / CLOCKS_PER_SEC;
	tmp = QString::number(cost);
	//out = tr("计算完成，共耗时") + tmp + "s，其中积分耗时" + QString::number(model->timer) + "s";
	out = tr("calculation complete.Total time:") + tmp + "s. Integration time:" + QString::number(model->timer) + "s";
	//out += "beta耗时" + QString::number(model->betaTimer) + "s，pre耗时" + QString::number(model->preTimer) + "s";
	//out += "子步" + QString::number(double(model->subSteps / model->endAndReversalPoint));
	//out += "。RK4的次数" + QString::number(double(model->CPM));
	message->close();
	delete message;
	ui.statusBar->showMessage(out);
}

void MainWindow::setFigureTitle() {
	*this->model->figureTitle = this->figureTitle->text();
}

void MainWindow::setAxisX(int x) {
	//ui.statusBar->showMessage(tr("设置X轴"));
	ui.statusBar->showMessage(tr("set axis X"));
	this->model->axisX = x;
}

void MainWindow::setAxisY(int y) {
	//ui.statusBar->showMessage(tr("设置Y轴"));
	ui.statusBar->showMessage(tr("set axis Y"));
	this->model->axisY = y;
}

void MainWindow::drawing() {
	//ui.statusBar->showMessage(tr("画图"));
	ui.statusBar->showMessage(tr("figure"));
	figure->setNewChart();
	figure->show();
}

void MainWindow::reset() {
	delete model, figure, displayParameter;
	model = new MODEL();
	figure = new Chart(model, this);
	displayParameter = new DisplayParameter(model, this);
	ui.statusBar->showMessage(tr("reset all parameter"));
}

void MainWindow::test() {
	QMessageBox *message = new QMessageBox(QMessageBox::Information, "测试中", tr("正在计算……"), QMessageBox::NoButton, this);
	message->show();
	filePath = tr("./para/test");
	resultFilePath = tr("./para/undrainedTriCompession");
	//resultFilePath = tr("./para/ansi");
	for (int i = 5; i < 8; i++) {
		filePath += QString::number(i) + tr(".dat");
		resultFilePath += QString::number(i) + tr(".txt");
		loadPara();
		model->Simulate();
		saveRes();
		reset();
		filePath = tr("./para/test");
		//resultFilePath = tr("./para/ansi");
		resultFilePath = tr("./para/undrainedTriCompession");
	}
	message->close();
}