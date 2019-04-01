#include "EndAndReversalStateNoCycle.h"

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

EndAndReversalStateNoCycle::EndAndReversalStateNoCycle(QWidget *parent) : QDialog(parent) {
	this->setWindowTitle("设置末状态");
	this->resize(330, 120);
	setFixedSize(this->width(), this->height());

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);

	labelType = new QLabel(this);
	labelType->setText("选择末状态类别");
	labelType->move(20, 20);
	labelType->resize(90, 30);

	labelPoint = new QLabel(this);
	labelPoint->setText("=");
	labelPoint->move(210, 20);
	labelPoint->resize(40, 30);

	inputType = new QComboBox(this);
	inputType->addItem("p");
	inputType->addItem("q");
	inputType->addItem("体应变");
	inputType->addItem("偏应变");
	inputType->addItem("步数");
	inputType->resize(80, 30);
	inputType->move(120, 20);

	inputPoint = new QLineEdit(this);
	inputPoint->resize(80, 30);
	inputPoint->move(230, 20);

	yes = new QPushButton(this);
	yes->setText("确定");
	yes->move(70, 75);
	yes->resize(60, 30);
	connect(yes, SIGNAL(clicked()), this, SLOT(clickYes()));

	cancel = new QPushButton(this);
	cancel->setText("取消");
	cancel->move(170, 75);
	cancel->resize(60, 30);
	connect(cancel, SIGNAL(clicked()), this, SLOT(clickCancel()));
}

EndAndReversalStateNoCycle::~EndAndReversalStateNoCycle() {
	delete labelPoint, labelType, inputPoint, inputType, yes, cancel;
}

void EndAndReversalStateNoCycle::clickYes() {
	int type, loop = 0;
	double point;
	type = inputType->currentIndex();
	point = inputPoint->text().toDouble();
	emit sendEndAndReversalState(type, point, loop);
	this->hide();
}

void EndAndReversalStateNoCycle::clickCancel() {
	this->hide();
}