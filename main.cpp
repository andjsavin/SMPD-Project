#include "mainwindow.h"
#include <QApplication>
#include <vectorcl.h>


int main(int argc, char *argv[])
{
    comb(64, 5);
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
