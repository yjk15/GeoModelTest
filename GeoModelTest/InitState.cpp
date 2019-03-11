#include "InitState.h"

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

InitState::InitState(QWidget *parent) : QDialog(parent) {
	this->resize(600, 230);
	this->setWindowTitle("设置初始应力应变");

	labelStress = new QLabel(this);
	labelStress->setText("初始应力(kPa)");
	labelStress->move(20, 20);
	labelStress->resize(90, 30);

	labelStrain = new QLabel(this);
	labelStrain->setText("初始应变");
	labelStrain->move(320, 20);
	labelStrain->resize(90, 30);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			inputStress[i * 3 + j].setParent(this);
			inputStress[i * 3 + j].setText("0");
			inputStress[i * 3 + j].resize(50, 30);
			inputStress[i * 3 + j].move(110 + 70 * j, 20 + 50 * i);
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			inputStrain[i * 3 + j].setParent(this);
			inputStrain[i * 3 + j].setText("0");
			inputStrain[i * 3 + j].resize(50, 30);
			inputStrain[i * 3 + j].move(390 + 70 * j, 20 + 50 * i);
		}
	}

	yes = new QPushButton(this);
	yes->setText("确定");
	yes->move(220, 180);
	yes->resize(60, 30);
	connect(yes, SIGNAL(clicked()), this, SLOT(clickYes()));

	cancel = new QPushButton(this);
	cancel->setText("取消");
	cancel->move(320, 180);
	cancel->resize(60, 30);
	connect(cancel, SIGNAL(clicked()), this, SLOT(clickCancel()));
}

InitState::~InitState() {
	delete yes, cancel, labelStress, labelStrain;
}

void InitState::clickYes() {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			stress(i, j) = inputStress[i * 3 + j].text().toDouble();
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			strain(i, j) = inputStrain[i * 3 + j].text().toDouble();
		}
	}

	emit sendInitState(stress, strain);
	this->hide();
}

void InitState::clickCancel() {
	this->hide();
}