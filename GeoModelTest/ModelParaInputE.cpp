#include "ModelParaInputE.h"

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

ModelParaInputE::ModelParaInputE(QWidget *parent) : QDialog(parent){
	this->resize(240, 190);
	this->setWindowTitle("设置线性模型参数");

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);

	labelE = new QLabel(this);
	labelE->setText("弹性模量E");
	labelE->move(20, 20);
	labelE->resize(80, 30);

	labelGPa = new QLabel(this);
	labelGPa->setText("kPa");
	labelGPa->move(200, 20);
	labelGPa->resize(80, 30);

	labelv = new QLabel(this);
	labelv->setText("泊松比E");
	labelv->move(20, 80);
	labelv->resize(80, 30);

	inputE = new QLineEdit(this);
	inputE->move(100, 20);
	inputE->resize(80, 30);

	inputv = new QLineEdit(this);
	inputv->move(100, 80);
	inputv->resize(80, 30);

	yes = new QPushButton(this);
	yes->setText("确定");
	yes->move(50, 140);
	yes->resize(60, 30);
	connect(yes, SIGNAL(clicked()), this, SLOT(clickYes()));

	cancel = new QPushButton(this);
	cancel->setText("取消");
	cancel->move(150, 140);
	cancel->resize(60, 30);
	connect(cancel, SIGNAL(clicked()), this, SLOT(closeDialog()));
}

ModelParaInputE::~ModelParaInputE() {
	delete labelE, labelGPa, labelv, E, v, yes, cancel;
}

void ModelParaInputE::clickYes() {
	E = inputE->text();
	v = inputv->text();
	emit sendEv(E.toDouble(), v.toDouble());
	this->hide();
}

void ModelParaInputE::closeDialog() {
	this->hide();
}