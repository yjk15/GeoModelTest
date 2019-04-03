#pragma once

#include <QtWidgets/QMainWindow>
#include <QString>
#include <ctime>
#include "ui_MainWindow.h"
#include "about.h"
#include "Model.h"
#include "chart.h"
#include "ModelParaInputE.h"
#include "ModelParaInputEB.h"
#include "ModelParaInputDM.h"
#include "ModelParaInputCycliq.h"
#include "InitState.h"
#include "EndAndReversalState.h"
#include "EndAndReversalStateNoCycle.h"
#include "DisplayParameter.h"

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = Q_NULLPTR);
	~MainWindow();

	void setMenu();
	void setType();
	void setModel();
	void calculate();
	void drawFigure();
	
	MODEL *model;
	ABOUT *about;
	QString filePath, resultFilePath;
	QLineEdit *figureTitle;
	Chart *figure;
	ModelParaInputE *inputModelE;
	ModelParaInputEB *inputModelEB;
	ModelParaInputDM *inputModelDM;
	ModelParaInputCycliq *inputModelCycliq;
	InitState *initState;
	EndAndReversalState *endAndReversalState;
	EndAndReversalStateNoCycle *endAndReversalStateNoCycle;
	DisplayParameter *displayParameter;

private:
	Ui::MainWindowClass ui;
	void loadPara();
	void saveRes();

private slots:
	void setInitState();
	void setEndState();
	void setTestType(int testType);
	void receiveInitState(MATRIX, MATRIX, double);
	void receiveEndAndReversalState(int, double, int);
	
	void openManual();
	void openAbout();
	void loadParameter();
	
	void saveParameter();
	void saveParameterInNewFile();
	void saveResult();
	
	void setConsModel(int model);
	void setModelParameter();
	void receiveE(double, double);
	void receiveEB(double, double, double, double, double);
	void receiveDM(double*);
	void receiveCycliq(double*);
	
	void setStepLength(int l);
	void checkParameter();
	void beginCalculate();

	void setFigureTitle();
	void setAxisX(int x);
	void setAxisY(int y);
	void drawing();
	void reset();

	void test();

};
