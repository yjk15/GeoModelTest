#include "ModelParaInputDMF.h"
#include <QtGui>

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

ModelParaInputDMF::ModelParaInputDMF(QWidget *parent) : QDialog(parent) {
	ui.setupUi(this);
	setFixedSize(this->width(), this->height());
	this->setWindowTitle("set DM with fabric parameter");
	input = new QLineEdit[16];

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);

	for (int i = 0; i < 16; i++) {
		input[i].setParent(this);
		input[i].move(260, 20 + 30 * i);
		input[i].resize(90, 21);
	}

	yes = new QPushButton(this);
	cancel = new QPushButton(this);

	yes->resize(75, 23);
	yes->move(80, 520);
	yes->setText("Yes");
	connect(yes, SIGNAL(clicked()), this, SLOT(clickYes()));

	cancel->resize(75, 23);
	cancel->move(230, 520);
	cancel->setText("Cancel");
	connect(cancel, SIGNAL(clicked()), this, SLOT(closeDialog()));
}

ModelParaInputDMF::~ModelParaInputDMF() {
	delete yes, cancel;
	delete[] input;
}

void ModelParaInputDMF::clickYes() {
	para[0] = 0;
	for (int i = 1; i < 17; i++) {
		para[i] = input[i - 1].text().toDouble();
	}

	emit sendParaDMF(para);
	this->hide();
}

void ModelParaInputDMF::closeDialog() {
	this->hide();
}