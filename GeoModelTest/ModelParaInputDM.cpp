#include "ModelParaInputDM.h"
#include <QtGui>

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

ModelParaInputDM::ModelParaInputDM(QWidget *parent) : QDialog(parent) {
	ui.setupUi(this);
	setFixedSize(this->width(), this->height());
	this->setWindowTitle("set DM parameter");
	input = new QLineEdit[15];
	method = new QComboBox(this);

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);

	for (int i = 0; i < 15; i++) {
		input[i].setParent(this);
		input[i].move(260, 20 + 30 * i);
		input[i].resize(90, 21);
	}
	method->resize(91, 21);
	method->move(160, 470);
	method->addItem("CPM");
	method->addItem("RK4");
	method->addItem("CPPM");

	yes = new QPushButton(this);
	cancel = new QPushButton(this);

	yes->resize(75, 23);
	yes->move(80, 510);
	yes->setText("Yes");
	connect(yes, SIGNAL(clicked()), this, SLOT(clickYes()));

	cancel->resize(75, 23);
	cancel->move(230, 510);
	cancel->setText("Cancel");
	connect(cancel, SIGNAL(clicked()), this, SLOT(closeDialog()));
}

ModelParaInputDM::~ModelParaInputDM() {
	delete yes, cancel, method;
	delete[] input;
}

void ModelParaInputDM::clickYes() {
	para[0] = method->currentIndex();
	for (int i = 1; i < 16; i++) {
		para[i] = input[i - 1].text().toDouble();
	}

	emit sendParaDM(para);
	this->hide();
}

void ModelParaInputDM::closeDialog() {
	this->hide();
}