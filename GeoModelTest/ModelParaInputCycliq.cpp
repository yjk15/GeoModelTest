#include "ModelParaInputCycliq.h"
#include <QtGui>

#pragma execution_character_set("utf-8")    // ��������������⣬ע�⣡����

ModelParaInputCycliq::ModelParaInputCycliq(QWidget *parent) : QDialog(parent) {
	ui.setupUi(this);
	setFixedSize(this->width(), this->height());
	this->setWindowTitle("����Cycliqģ�Ͳ���");
	input = new QLineEdit[14];
	//method = new QComboBox(this);

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);

	for (int i = 0; i < 14; i++) {
		input[i].setParent(this);
		input[i].move(120, 20 + 30 * i);
		input[i].resize(90, 21);
	}
	//method->resize(91, 21);
	//method->move(160, 470);
	//method->addItem("��ʽ����");
	//method->addItem("��ʽ����");

	yes = new QPushButton(this);
	cancel = new QPushButton(this);

	yes->resize(75, 23);
	yes->move(40, 450);
	yes->setText("ȷ��");
	connect(yes, SIGNAL(clicked()), this, SLOT(clickYes()));

	cancel->resize(75, 23);
	cancel->move(140, 450);
	cancel->setText("ȡ��");
	connect(cancel, SIGNAL(clicked()), this, SLOT(closeDialog()));
}

ModelParaInputCycliq::~ModelParaInputCycliq() {
	delete yes, cancel;// method;
	delete[] input;
}

void ModelParaInputCycliq::clickYes() {
	para[0] = 0;
	for (int i = 1; i < 15; i++) {
		para[i] = input[i - 1].text().toDouble();
	}

	emit sendParaCycliq(para);
	this->hide();
}

void ModelParaInputCycliq::closeDialog() {
	this->hide();
}