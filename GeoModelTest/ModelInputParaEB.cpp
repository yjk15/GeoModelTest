#include "ModelParaInputEB.h"

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

ModelParaInputEB::ModelParaInputEB(QWidget *parent) : QDialog(parent) {
	this->resize(240, 330);
	setFixedSize(this->width(), this->height());
	this->setWindowTitle("设置EB模型参数");

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);

	label0 = new QLabel(this);
	label0->resize(60, 30);
	label0->move(20, 20);
	label0->setText("Kd1");

	label1 = new QLabel(this);
	label1->resize(60, 30);
	label1->move(20, 70);
	label1->setText("Kd2");

	label2 = new QLabel(this);
	label2->resize(60, 30);
	label2->move(20, 120);
	label2->setText("nd2");

	label3 = new QLabel(this);
	label3->resize(60, 30);
	label3->move(20, 170);
	label3->setText("nud");

	label4 = new QLabel(this);
	label4->resize(60, 30);
	label4->move(20, 220);
	label4->setText("gammamax");
	
	inputEB0 = new QLineEdit(this);
	inputEB0->move(130, 20);
	inputEB0->resize(90, 30);

	inputEB1 = new QLineEdit(this);
	inputEB1->move(130, 70);
	inputEB1->resize(90, 30);

	inputEB2 = new QLineEdit(this);
	inputEB2->move(130, 120);
	inputEB2->resize(90, 30);

	inputEB3 = new QLineEdit(this);
	inputEB3->move(130, 170);
	inputEB3->resize(90, 30);

	inputEB4 = new QLineEdit(this);
	inputEB4->move(130, 220);
	inputEB4->resize(90, 30);

	yes = new QPushButton(this);
	yes->setText("确定");
	yes->move(50, 270);
	yes->resize(60, 30);
	connect(yes, SIGNAL(clicked()), this, SLOT(clickYes()));

	cancel = new QPushButton(this);
	cancel->setText("取消");
	cancel->move(150, 270);
	cancel->resize(60, 30);
	connect(cancel, SIGNAL(clicked()), this, SLOT(closeDialog()));
}

ModelParaInputEB::~ModelParaInputEB() {
	delete label0, label1, label2, label3, label4, inputEB0, inputEB1, inputEB2, inputEB3, inputEB4;
	delete yes, cancel;
}

void ModelParaInputEB::clickYes() {
	Kd1 = inputEB0->text();
	Kd2 = inputEB1->text();
	nd2 = inputEB2->text();
	nud = inputEB3->text();
	gammamax = inputEB4->text();
	emit sendEBPara(Kd1.toDouble(), Kd2.toDouble(), nd2.toDouble(), nud.toDouble(), gammamax.toDouble());
	this->hide();
}

void ModelParaInputEB::closeDialog() {
	this->hide();
}