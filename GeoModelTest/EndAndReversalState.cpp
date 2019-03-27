#include "EndAndReversalState.h"

#pragma execution_character_set("utf-8")    // ��������������⣬ע�⣡����

EndAndReversalState::EndAndReversalState(QWidget *parent) : QDialog(parent) {
	this->setWindowTitle("����ĩ״̬");
	this->resize(330, 200);

	labelType = new QLabel(this);
	labelType->setText("ѡ��ĩ״̬���");
	labelType->move(20, 20);
	labelType->resize(90, 30);

	labelPoint = new QLabel(this);
	labelPoint->setText("=");
	labelPoint->move(210, 20);
	labelPoint->resize(40, 30);

	inputType = new QComboBox(this);
	inputType->addItem("p");
	inputType->addItem("q");
	inputType->addItem("��Ӧ��");
	inputType->addItem("ƫӦ��");
	inputType->resize(80, 30);
	inputType->move(120, 20);

	inputPoint = new QLineEdit(this);
	inputPoint->resize(80, 30);
	inputPoint->move(230, 20);

	labelLoop = new QLabel(this);
	labelLoop->setText("��ת����");
	labelLoop->move(20, 90);
	labelLoop->resize(60, 30);

	inputLoop = new QLineEdit(this);
	inputLoop->setText("0");
	inputLoop->resize(80, 30);
	inputLoop->move(120, 90);

	yes = new QPushButton(this);
	yes->setText("ȷ��");
	yes->move(70, 150);
	yes->resize(60, 30);
	connect(yes, SIGNAL(clicked()), this, SLOT(clickYes()));

	cancel = new QPushButton(this);
	cancel->setText("ȡ��");
	cancel->move(170, 150);
	cancel->resize(60, 30);
	connect(cancel, SIGNAL(clicked()), this, SLOT(clickCancel()));
}

EndAndReversalState::~EndAndReversalState() {
	delete labelPoint, labelType, inputPoint, inputType, yes, cancel;
}

void EndAndReversalState::clickYes() {
	int type, loop;
	double point;
	type = inputType->currentIndex();
	point = inputPoint->text().toDouble();
	loop = inputLoop->text().toInt();
	emit sendEndAndReversalState(type, point, loop);
	this->hide();
}

void EndAndReversalState::clickCancel() {
	this->hide();
}