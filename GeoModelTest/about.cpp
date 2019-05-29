#include "about.h"

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

ABOUT::ABOUT(QWidget *parent) : QDialog(parent){
	//初始化对话框的大小以及设置其名字
	this->resize(320, 110);
	setFixedSize(this->width(), this->height());
	this->setWindowTitle("about");

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);


	label1 = new QLabel(this);
	label2 = new QLabel(this);
	label3 = new QLabel(this);
	label1->setText("This software is only used for study.");
	label1->move(20, 20);
	label2->setText("Github address is https://github.com/yjk15/GeoModelTest");
	label2->move(20, 50);
	label3->setText("Thanks for the help for the teacher.");
	label3->move(20, 80);
}

ABOUT::~ABOUT() {}