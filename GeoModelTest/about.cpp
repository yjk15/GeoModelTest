#include "about.h"

#pragma execution_character_set("utf-8")    // ��������������⣬ע�⣡����

ABOUT::ABOUT(QWidget *parent) : QDialog(parent){
	//��ʼ���Ի���Ĵ�С�Լ�����������
	this->resize(280, 110);
	this->setWindowTitle("����");

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);


	label1 = new QLabel(this);
	label2 = new QLabel(this);
	label3 = new QLabel(this);
	label1->setText("�����������ѧϰ����������������ҵ��;��");
	label1->move(20, 20);
	label2->setText("github��ַΪ");
	label2->move(20, 50);
	label3->setText("��л��ʦ��ͬѧ�İ�����");
	label3->move(20, 80);
}

ABOUT::~ABOUT() {}