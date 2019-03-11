/********************************************************************************
** Form generated from reading UI file 'chart.ui'
**
** Created by: Qt User Interface Compiler version 5.13.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CHART_H
#define UI_CHART_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_chartClass
{
public:
    QWidget *centralWidget;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *chartClass)
    {
        if (chartClass->objectName().isEmpty())
            chartClass->setObjectName(QString::fromUtf8("chartClass"));
        chartClass->resize(600, 400);
        centralWidget = new QWidget(chartClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        chartClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(chartClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 600, 23));
        chartClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(chartClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        chartClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(chartClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        chartClass->setStatusBar(statusBar);

        retranslateUi(chartClass);

        QMetaObject::connectSlotsByName(chartClass);
    } // setupUi

    void retranslateUi(QMainWindow *chartClass)
    {
        chartClass->setWindowTitle(QCoreApplication::translate("chartClass", "chart", nullptr));
    } // retranslateUi

};

namespace Ui {
    class chartClass: public Ui_chartClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CHART_H
