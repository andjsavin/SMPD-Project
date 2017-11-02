#include "mainwindow.h"
#include <QApplication>
#include <vectorcl.h>


int main(int argc, char *argv[])
{
    std::vector<std::vector<int>> v = comb(5, 3);
    printv(v);
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
