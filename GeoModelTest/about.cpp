#include "about.h"

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

ABOUT::ABOUT(QWidget *parent) : QDialog(parent){
	//初始化对话框的大小以及设置其名字
	this->resize(280, 110);
	this->setWindowTitle("关于");

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);


	label1 = new QLabel(this);
	label2 = new QLabel(this);
	label3 = new QLabel(this);
	label1->setText("本软件仅用于学习、交流，不用做商业用途。");
	label1->move(20, 20);
	label2->setText("github地址为");
	label2->move(20, 50);
	label3->setText("感谢老师、同学的帮助。");
	label3->move(20, 80);
}

ABOUT::~ABOUT() {}