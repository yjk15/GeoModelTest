#include "chart.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	chart w;
	w.show();
	return a.exec();
}
