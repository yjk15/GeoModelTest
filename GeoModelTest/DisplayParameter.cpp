#include "DisplayParameter.h"

#pragma execution_character_set("utf-8")    // ��������������⣬ע�⣡����

DisplayParameter::DisplayParameter(MODEL *m, QWidget *parent) : QDialog(parent) {
	//��ʼ���Ի���Ĵ�С�Լ�����������
	this->resize(600, 600);
	this->setWindowTitle("������");

	model = m;
	label = new QLabel[50];
	DisplayPara(model->model, model->testType);
}

DisplayParameter::~DisplayParameter() {
	delete[] label;
}

void DisplayParameter::DisplayPara(int modelType, int testType) {
	label[0].setParent(this);
	label[0].setText("�������");
	label[0].move(20, 20);
	label[0].resize(90, 30);

	label[1].setParent(this);
	switch (model->testType){
	case 0:
		label[1].setText("����ˮ����ѹ������");
		break;
	case 1:
		label[1].setText("����ˮ���ἷ������");
		break;
	case 2:
		label[1].setText("����ˮ����ѭ������");
		break;
	case 3:
		label[1].setText("��ˮ����ѹ������");
		break;
	case 4:
		label[1].setText("��ˮ���ἷ������");
		break;
	case 5:
		label[1].setText("��ˮ����ѭ������");
		break;
	default:
		break;
	}
	label[1].move(80, 20);
	label[1].resize(200, 30);

	label[2].setParent(this);
	label[2].setText("��ʼӦ��(kPa)");
	label[2].move(20, 120);
	label[2].resize(90, 30);

	label[3].setParent(this);
	label[3].setText("��ʼӦ��");
	label[3].move(320, 120);
	label[3].resize(90, 30);

	QString tmp;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			label[i * 3 + j + 4].setParent(this);
			tmp = QString::number(model->stress(i, j));
			label[i * 3 + j + 4].setText(tmp);
			label[i * 3 + j + 4].resize(50, 30);
			label[i * 3 + j + 4].move(110 + 70 * j, 70 + 50 * i);
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			label[i * 3 + j + 13].setParent(this);
			tmp = QString::number(model->strain(i, j));
			label[i * 3 + j + 13].setText(tmp);
			label[i * 3 + j + 13].resize(50, 30);
			label[i * 3 + j + 13].move(390 + 70 * j, 70 + 50 * i);
		}
	}

	label[22].setParent(this);
	label[22].setText("ĩ״̬���");
	label[22].move(20, 220);
	label[22].resize(90, 30);

	label[23].setParent(this);
	label[23].setText("=");
	label[23].move(180, 220);
	label[23].resize(40, 30);

	label[24].setParent(this);
	switch (model->endAndReversalType) {
	case 0:
		label[24].setText("p");
		break;
	case 1:
		label[24].setText("q");
		break;
	case 2:
		label[24].setText("��Ӧ��");
		break;
	default:
		break;
	}
	label[24].resize(80, 30);
	label[24].move(120, 220);

	label[25].setParent(this);
	tmp = QString::number(model->endAndReversalPoint);
	label[25].setText(tmp);
	label[25].resize(80, 30);
	label[25].move(230, 220);

	if (model->testType == 2 || model->testType == 5) {
		label[26].setParent(this);
		label[26].setText("ѭ������");
		label[26].move(20, 270);
		label[26].resize(60, 30);

		label[27].setParent(this);
		label[27].setText("0");
		label[27].resize(80, 30);
		label[27].move(120, 270);
	}

	label[28].setParent(this);
	label[28].setText("����ģ��");
	label[28].resize(80, 30);
	label[28].move(20, 320);

	label[29].setParent(this);
	switch (model->model) {
	case 0:
		label[29].setText("����ģ��");
		break;
	case 1:
		label[29].setText("EBģ��");
		break;
	case 2:
		label[29].setText("Dafalias and Manzariģ��");
		break;
	default:
		break;
	}
	label[29].resize(200, 30);
	label[29].move(80, 320);
	switch (modelType) {
	case 0:
		label[30].setParent(this);
		label[30].setText("����ģ��");
		label[30].move(20, 370);
		label[30].resize(90, 30);

		label[31].setParent(this);
		tmp = QString::number(model->internalParameter[0]);
		label[31].setText(tmp);
		label[31].move(80, 370);
		label[31].resize(90, 30);

		label[32].setParent(this);
		label[32].setText("���ɱ�");
		label[32].move(20, 420);
		label[32].resize(90, 30);

		label[33].setParent(this);
		tmp = QString::number(model->internalParameter[1]);
		label[33].setText(tmp);
		label[33].move(80, 420);
		label[33].resize(90, 30);
		break;
	default:
		break;
	}
}