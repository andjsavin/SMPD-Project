#include "mainwindow.h"
#include <QApplication>
#include <vectorcl.h>


int main(int argc, char *argv[])
{
    comb(10, 4);
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
