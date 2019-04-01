#include "DisplayParameter.h"

#pragma execution_character_set("utf-8")    // 解决汉字乱码问题，注意！！！

DisplayParameter::DisplayParameter(MODEL *m, QWidget *parent) : QDialog(parent) {
	//初始化对话框的大小以及设置其名字
	this->resize(600, 660);
	setFixedSize(this->width(), this->height());
	this->setWindowTitle("检查参数");

	Qt::WindowFlags flags = Qt::Dialog;
	flags |= Qt::WindowCloseButtonHint;
	setWindowFlags(flags);

	model = m;
	label = new QLabel[70];
	DisplayPara(model->model, model->testType);
}

DisplayParameter::~DisplayParameter() {
	delete[] label;
}

void DisplayParameter::DisplayPara(int modelType, int testType) {
	QString tmp;

	label[0].setParent(this);
	label[0].setText("试验类别");
	label[0].move(20, 20);
	label[0].resize(90, 30);

	label[1].setParent(this);
	switch (model->testType){
	case 0:
		label[1].setText("不排水三轴压缩试验");
		break;
	case 1:
		label[1].setText("不排水三轴挤长试验");
		break;
	case 2:
		label[1].setText("不排水三轴循环试验");
		break;
	case 3:
		label[1].setText("排水三轴压缩试验");
		break;
	case 4:
		label[1].setText("排水三轴挤长试验");
		break;
	case 5:
		label[1].setText("排水三轴循环试验");
		break;
	default:
		break;
	}
	label[1].move(80, 20);
	label[1].resize(200, 30);

	labelee[0].setParent(this);
	labelee[0].setText("模拟开始时孔隙比e");
	labelee[0].move(300, 20);
	labelee[0].resize(120, 30);

	labelee[1].setParent(this);
	tmp = QString::number(model->ee);
	labelee[1].setText(tmp);
	labelee[1].move(420, 20);
	labelee[1].resize(90, 30);

	labelee[2].setParent(this);
	labelee[2].setText("计算步长");
	labelee[2].move(300, 320);
	labelee[2].resize(90, 30);

	labelee[3].setParent(this);
	tmp = QString::number(model->stepLength);
	labelee[3].setText(tmp);
	labelee[3].move(350, 320);
	labelee[3].resize(90, 30);

	label[2].setParent(this);
	label[2].setText("初始应力(kPa)");
	label[2].move(20, 120);
	label[2].resize(90, 30);

	label[3].setParent(this);
	label[3].setText("初始应变");
	label[3].move(320, 120);
	label[3].resize(90, 30);
	
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
	label[22].setText("末状态类别");
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
		label[24].setText("体应变");
		break;
	case 3:
		label[24].setText("偏应变");
		break;
	case 4:
		label[24].setText("步数");
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
		label[26].setText("反转次数");
		label[26].move(20, 270);
		label[26].resize(60, 30);

		label[27].setParent(this);
		tmp = QString::number(model->loop);
		label[27].setText(tmp);
		label[27].resize(80, 30);
		label[27].move(120, 270);
	}

	label[28].setParent(this);
	label[28].setText("本构模型");
	label[28].resize(80, 30);
	label[28].move(20, 320);

	label[29].setParent(this);
	switch (model->model) {
	case 0:
		label[29].setText("线性模型");
		break;
	case 1:
		label[29].setText("EB模型");
		break;
	case 2:
		label[29].setText("Dafalias and Manzari模型");
		break;
	default:
		break;
	}
	label[29].resize(200, 30);
	label[29].move(80, 320);
	switch (modelType) {
	case 0:
		label[30].setParent(this);
		label[30].setText("弹性模量");
		label[30].move(20, 370);
		label[30].resize(90, 30);

		label[31].setParent(this);
		tmp = QString::number(model->internalParameter[0]);
		label[31].setText(tmp);
		label[31].move(80, 370);
		label[31].resize(90, 30);

		label[32].setParent(this);
		label[32].setText("泊松比");
		label[32].move(20, 420);
		label[32].resize(90, 30);

		label[33].setParent(this);
		tmp = QString::number(model->internalParameter[1]);
		label[33].setText(tmp);
		label[33].move(80, 420);
		label[33].resize(90, 30);
		break;
	case 1:
		label[30].setParent(this);
		label[30].setText("Kd1");
		label[30].move(20, 370);
		label[30].resize(90, 30);

		label[31].setParent(this);
		tmp = QString::number(model->internalParameter[0]);
		label[31].setText(tmp);
		label[31].move(80, 370);
		label[31].resize(90, 30);

		label[32].setParent(this);
		label[32].setText("Kd2");
		label[32].move(20, 420);
		label[32].resize(90, 30);

		label[33].setParent(this);
		tmp = QString::number(model->internalParameter[1]);
		label[33].setText(tmp);
		label[33].move(80, 420);
		label[33].resize(90, 30);

		label[34].setParent(this);
		label[34].setText("nd2");
		label[34].move(20, 470);
		label[34].resize(90, 30);

		label[35].setParent(this);
		tmp = QString::number(model->internalParameter[2]);
		label[35].setText(tmp);
		label[35].move(80, 470);
		label[35].resize(90, 30);

		label[36].setParent(this);
		label[36].setText("nud");
		label[36].move(20, 520);
		label[36].resize(90, 30);

		label[37].setParent(this);
		tmp = QString::number(model->internalParameter[3]);
		label[37].setText(tmp);
		label[37].move(80, 520);
		label[37].resize(90, 30);

		label[38].setParent(this);
		label[38].setText("gammamax");
		label[38].move(20, 570);
		label[38].resize(90, 30);

		label[39].setParent(this);
		tmp = QString::number(model->internalParameter[4]);
		label[39].setText(tmp);
		label[39].move(80, 570);
		label[39].resize(90, 30);
		break;
	case 2:
		label[30].setParent(this);
		label[30].setText("G_0");
		label[30].move(20, 370);
		label[30].resize(90, 30);

		label[31].setParent(this);
		tmp = QString::number(model->internalParameter[1]);
		label[31].setText(tmp);
		label[31].move(80, 370);
		label[31].resize(90, 30);

		label[32].setParent(this);
		label[32].setText("v");
		label[32].move(250, 370);
		label[32].resize(90, 30);

		label[33].setParent(this);
		tmp = QString::number(model->internalParameter[2]);
		label[33].setText(tmp);
		label[33].move(310, 370);
		label[33].resize(90, 30);

		label[34].setParent(this);
		label[34].setText("M");
		label[34].move(460, 370);
		label[34].resize(90, 30);

		label[35].setParent(this);
		tmp = QString::number(model->internalParameter[3]);
		label[35].setText(tmp);
		label[35].move(520, 370);
		label[35].resize(90, 30);

		label[36].setParent(this);
		label[36].setText("c");
		label[36].move(20, 420);
		label[36].resize(90, 30);

		label[37].setParent(this);
		tmp = QString::number(model->internalParameter[4]);
		label[37].setText(tmp);
		label[37].move(80, 420);
		label[37].resize(90, 30);

		label[38].setParent(this);
		label[38].setText("λc");
		label[38].move(250, 420);
		label[38].resize(90, 30);

		label[39].setParent(this);
		tmp = QString::number(model->internalParameter[5]);
		label[39].setText(tmp);
		label[39].move(310, 420);
		label[39].resize(90, 30);

		label[40].setParent(this);
		label[40].setText("e_0");
		label[40].move(460, 420);
		label[40].resize(90, 30);

		label[41].setParent(this);
		tmp = QString::number(model->internalParameter[6]);
		label[41].setText(tmp);
		label[41].move(520, 420);
		label[41].resize(90, 30);

		label[42].setParent(this);
		label[42].setText("ξ");
		label[42].move(20, 470);
		label[42].resize(90, 30);

		label[43].setParent(this);
		tmp = QString::number(model->internalParameter[7]);
		label[43].setText(tmp);
		label[43].move(80, 470);
		label[43].resize(90, 30);

		label[44].setParent(this);
		label[44].setText("m");
		label[44].move(250, 470);
		label[44].resize(90, 30);

		label[45].setParent(this);
		tmp = QString::number(model->internalParameter[8]);
		label[45].setText(tmp);
		label[45].move(310, 470);
		label[45].resize(90, 30);

		label[46].setParent(this);
		label[46].setText("h_0");
		label[46].move(460, 470);
		label[46].resize(90, 30);

		label[47].setParent(this);
		tmp = QString::number(model->internalParameter[9]);
		label[47].setText(tmp);
		label[47].move(520, 470);
		label[47].resize(90, 30);

		label[48].setParent(this);
		label[48].setText("c_h");
		label[48].move(20, 520);
		label[48].resize(90, 30);

		label[49].setParent(this);
		tmp = QString::number(model->internalParameter[10]);
		label[49].setText(tmp);
		label[49].move(80, 520);
		label[49].resize(90, 30);

		label[50].setParent(this);
		label[50].setText("n^b");
		label[50].move(250, 520);
		label[50].resize(90, 30);

		label[51].setParent(this);
		tmp = QString::number(model->internalParameter[11]);
		label[51].setText(tmp);
		label[51].move(310, 520);
		label[51].resize(90, 30);

		label[52].setParent(this);
		label[52].setText("A_0");
		label[52].move(460, 520);
		label[52].resize(90, 30);

		label[53].setParent(this);
		tmp = QString::number(model->internalParameter[12]);
		label[53].setText(tmp);
		label[53].move(520, 520);
		label[53].resize(90, 30);

		label[54].setParent(this);
		label[54].setText("n^d");
		label[54].move(20, 570);
		label[54].resize(90, 30);

		label[55].setParent(this);
		tmp = QString::number(model->internalParameter[13]);
		label[55].setText(tmp);
		label[55].move(80, 570);
		label[55].resize(90, 30);

		label[56].setParent(this);
		label[56].setText("z_{max}");
		label[56].move(250, 570);
		label[56].resize(90, 30);

		label[57].setParent(this);
		tmp = QString::number(model->internalParameter[14]);
		label[57].setText(tmp);
		label[57].move(310, 570);
		label[57].resize(90, 30);

		label[58].setParent(this);
		label[58].setText("c_z");
		label[58].move(460, 570);
		label[58].resize(90, 30);

		label[59].setParent(this);
		tmp = QString::number(model->internalParameter[15]);
		label[59].setText(tmp);
		label[59].move(520, 570);
		label[59].resize(90, 30);

		label[60].setParent(this);
		label[60].setText("积分方法");
		label[60].move(20, 620);
		label[60].resize(90, 30);

		label[61].setParent(this);
		if (model->internalParameter[0] == 0)
			label[61].setText("隐式积分");
		if (model->internalParameter[0] == 1)
			label[61].setText("显式积分");
		label[61].move(80, 620);
		label[61].resize(90, 30);
		break;
	case 3:
		label[30].setParent(this);
		label[30].setText("G_0");
		label[30].move(20, 370);
		label[30].resize(90, 30);

		label[31].setParent(this);
		tmp = QString::number(model->internalParameter[1]);
		label[31].setText(tmp);
		label[31].move(80, 370);
		label[31].resize(90, 30);

		label[32].setParent(this);
		label[32].setText("κ");
		label[32].move(250, 370);
		label[32].resize(90, 30);

		label[33].setParent(this);
		tmp = QString::number(model->internalParameter[2]);
		label[33].setText(tmp);
		label[33].move(310, 370);
		label[33].resize(90, 30);

		label[34].setParent(this);
		label[34].setText("h");
		label[34].move(460, 370);
		label[34].resize(90, 30);

		label[35].setParent(this);
		tmp = QString::number(model->internalParameter[3]);
		label[35].setText(tmp);
		label[35].move(520, 370);
		label[35].resize(90, 30);

		label[36].setParent(this);
		label[36].setText("M");
		label[36].move(20, 420);
		label[36].resize(90, 30);

		label[37].setParent(this);
		tmp = QString::number(model->internalParameter[4]);
		label[37].setText(tmp);
		label[37].move(80, 420);
		label[37].resize(90, 30);

		label[38].setParent(this);
		label[38].setText("d_{re,1}");
		label[38].move(250, 420);
		label[38].resize(90, 30);

		label[39].setParent(this);
		tmp = QString::number(model->internalParameter[5]);
		label[39].setText(tmp);
		label[39].move(310, 420);
		label[39].resize(90, 30);

		label[40].setParent(this);
		label[40].setText("d_{re,2}");
		label[40].move(460, 420);
		label[40].resize(90, 30);

		label[41].setParent(this);
		tmp = QString::number(model->internalParameter[6]);
		label[41].setText(tmp);
		label[41].move(520, 420);
		label[41].resize(90, 30);

		label[42].setParent(this);
		label[42].setText("d_{ir}");
		label[42].move(20, 470);
		label[42].resize(90, 30);

		label[43].setParent(this);
		tmp = QString::number(model->internalParameter[7]);
		label[43].setText(tmp);
		label[43].move(80, 470);
		label[43].resize(90, 30);

		label[44].setParent(this);
		label[44].setText("α");
		label[44].move(250, 470);
		label[44].resize(90, 30);

		label[45].setParent(this);
		tmp = QString::number(model->internalParameter[8]);
		label[45].setText(tmp);
		label[45].move(310, 470);
		label[45].resize(90, 30);

		label[46].setParent(this);
		label[46].setText("γ_{d,r}");
		label[46].move(460, 470);
		label[46].resize(90, 30);

		label[47].setParent(this);
		tmp = QString::number(model->internalParameter[9]);
		label[47].setText(tmp);
		label[47].move(520, 470);
		label[47].resize(90, 30);

		label[48].setParent(this);
		label[48].setText("n^p");
		label[48].move(20, 520);
		label[48].resize(90, 30);

		label[49].setParent(this);
		tmp = QString::number(model->internalParameter[10]);
		label[49].setText(tmp);
		label[49].move(80, 520);
		label[49].resize(90, 30);

		label[50].setParent(this);
		label[50].setText("n^d");
		label[50].move(250, 520);
		label[50].resize(90, 30);

		label[51].setParent(this);
		tmp = QString::number(model->internalParameter[11]);
		label[51].setText(tmp);
		label[51].move(310, 520);
		label[51].resize(90, 30);

		label[52].setParent(this);
		label[52].setText("λ_c");
		label[52].move(460, 520);
		label[52].resize(90, 30);

		label[53].setParent(this);
		tmp = QString::number(model->internalParameter[12]);
		label[53].setText(tmp);
		label[53].move(520, 520);
		label[53].resize(90, 30);

		label[54].setParent(this);
		label[54].setText("e_0");
		label[54].move(20, 570);
		label[54].resize(90, 30);

		label[55].setParent(this);
		tmp = QString::number(model->internalParameter[13]);
		label[55].setText(tmp);
		label[55].move(80, 570);
		label[55].resize(90, 30);

		label[56].setParent(this);
		label[56].setText("ξ");
		label[56].move(250, 570);
		label[56].resize(90, 30);

		label[57].setParent(this);
		tmp = QString::number(model->internalParameter[14]);
		label[57].setText(tmp);
		label[57].move(310, 570);
		label[57].resize(90, 30);
		break;
	default:
		break;
	}
}