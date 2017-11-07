#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include "vectorcl.h"
#include "matrixutil.hpp"


#include <QImage>
#include <QDebug>




MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    FSupdateButtonState();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::updateDatabaseInfo()
{
    ui->FScomboBox->clear();
    for(unsigned int i=1; i<=database.getNoFeatures(); ++i)
        ui->FScomboBox->addItem(QString::number(i));

    ui->FStextBrowserDatabaseInfo->setText("noClass: " +  QString::number(database.getNoClass()));
    ui->FStextBrowserDatabaseInfo->append("noObjects: "  +  QString::number(database.getNoObjects()));
    ui->FStextBrowserDatabaseInfo->append("noFeatures: "  +  QString::number(database.getNoFeatures()));

}

void MainWindow::FSupdateButtonState(void)
{
    if(database.getNoObjects()==0)
    {
        FSsetButtonState(false);
    }
    else
        FSsetButtonState(true);

}


void MainWindow::FSsetButtonState(bool state)
{
   ui->FScomboBox->setEnabled(state);
   ui->FSpushButtonCompute->setEnabled(state);
   ui->FSpushButtonSaveFile->setEnabled(state);
   ui->FSradioButtonFisher->setEnabled(state);
   ui->FSradioButtonSFS->setEnabled(state);
}

void MainWindow::on_FSpushButtonOpenFile_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Open TextFile"), "", tr("Texts Files (*.txt)"));

    if ( !database.load(fileName.toStdString()) )
        QMessageBox::warning(this, "Warning", "File corrupted !!!");
    else
        QMessageBox::information(this, fileName, "File loaded !!!");

    FSupdateButtonState();
    updateDatabaseInfo();
}

void MainWindow::on_FSpushButtonCompute_clicked()
{
    int dimension = ui->FScomboBox->currentText().toInt();


    if( ui->FSradioButtonFisher ->isChecked())
    {
        if (database.getNoClass() == 2)
        {
            std::vector<Object> class1;
            std::vector<Object> class2;
            std::vector<Object> all_obj = database.getObjects();
            std::string cur_name = database.getObjects()[0].getClassName();
            for (uint i = 0; i < all_obj.size(); i++) {
                if (all_obj[i].getClassName() == cur_name) {
                    class1.push_back(all_obj[i]);
                } else {
                    class2.push_back(all_obj[i]);
                }
            }
            std::vector<std::vector<double>> class1ObFeatures;
            std::vector<std::vector<double>> class2ObFeatures;
            for (uint i = 0; i < class1.size(); i++) {
                class1ObFeatures.push_back(class1[i].getFeatures());
            }
            for (uint i = 0; i < class2.size(); i++) {
                class2ObFeatures.push_back(class2[i].getFeatures());
            }
            class1ObFeatures = transponate(class1ObFeatures);
            class2ObFeatures = transponate(class2ObFeatures);
            std::vector<std::vector<double>> cl1MAvg = minusAvg(getMatrixMedian(class1ObFeatures), class1ObFeatures);
            std::vector<std::vector<double>> cl2MAvg = minusAvg(getMatrixMedian(class2ObFeatures), class2ObFeatures);
            int N = database.getNoFeatures();
            int K = dimension;
            std::string bitmask(K, 1);
            bitmask.resize(N, 0);
            double max = 0.0;
            std::vector<int> maxVector;
            do {
                std::vector<int> vh;
                for (int i = 0; i < N; ++i)
                {
                    if (bitmask[i]) {
                        vh.push_back(i);
                    }
                }
                std::vector<std::vector<double>> comb1 = getMatrixFromVector(vh, cl1MAvg);
                std::vector<std::vector<double>> comb2 = getMatrixFromVector(vh, cl2MAvg);
                std::vector<std::vector<double>> mm1 = multiplyMatrix(comb1, transponate(comb1), getProbabilityVector(class1.size()));
                std::vector<std::vector<double>> mm2 = multiplyMatrix(comb2, transponate(comb2), getProbabilityVector(class2.size()));
                std::vector<double> median1 = getVectorFromVector(vh, getMatrixMedian(class1ObFeatures));
                std::vector<double> median2 = getVectorFromVector(vh, getMatrixMedian(class2ObFeatures));
                double medianModule = getVectorModule(getVectorDifference(median1, median2));
                bnu::matrix<double> m1 = getMatrix(mm1);
                bnu::matrix<double> m2 = getMatrix(mm2);
                double dt1 = determinant(m1);
                double dt2 = determinant(m2);
                if (dimension == 1) {
                    if (max < medianModule/(sqrt(dt1) + sqrt(dt2)))
                    {
                        maxVector = vh;
                        max = medianModule/(sqrt(dt1) + sqrt(dt2));
                    }
                } else {
                    if (max < medianModule/(dt1 + dt2))
                    {
                        maxVector= vh;
                        max = medianModule/(dt1 + dt2);
                    }
                }
            } while (std::prev_permutation(bitmask.begin(), bitmask.end())); // lexicographilly permute bitmask
            QString s = vectorToString(maxVector);
            ui->FStextBrowserDatabaseInfo->append("max_vector: {"  +  s + "} " + QString::number(max));
//            std::vector<std::vector<int>> featureCombinations = comb(database.getNoFeatures(), dimension);
//            std::map<std::vector<int>, double> ftrs;
//            for (int i = 0; i < featureCombinations.size(); i++) {
//                std::vector<std::vector<double>> comb1 = getMatrixFromVector(featureCombinations[i], cl1MAvg);
//                std::vector<std::vector<double>> comb2 = getMatrixFromVector(featureCombinations[i], cl2MAvg);
//                std::vector<std::vector<double>> mm1 = multiplyMatrix(comb1, transponate(comb1), getProbabilityVector(class1.size()));
//                std::vector<std::vector<double>> mm2 = multiplyMatrix(comb2, transponate(comb2), getProbabilityVector(class2.size()));
//                std::vector<double> median1 = getVectorFromVector(featureCombinations[i], getMatrixMedian(class1ObFeatures));
//                std::vector<double> median2 = getVectorFromVector(featureCombinations[i], getMatrixMedian(class2ObFeatures));
//                double medianModule = getVectorModule(getVectorDifference(median1, median2));
//                double dt1 = getMatrixDeterminant(mm1);
//                double dt2 = getMatrixDeterminant(mm2);
//                if (dimension == 1) {
//                    ftrs[featureCombinations[i]] = medianModule/(sqrt(dt1) + sqrt(dt2));
//                } else {
//                ftrs[featureCombinations[i]] = medianModule/(dt1 + dt2);
//                }
//            }
//            double max = 0.0;
//            std::pair<std::vector<int>, double> maxp;
//            for (auto it = ftrs.begin(); it != ftrs.end(); ++it) {
//                if (it->second > max) {
//                    max = it->second;
//                    maxp = std::make_pair(it->first, it->second);
//                }
//            }
//            QString s = vectorToString(maxp.first);
//            ui->FStextBrowserDatabaseInfo->append("max_vector: {"  +  s + "} " + QString::number(maxp.second));
//        {
//            double FLD = 0, tmp;
//            int max_ind = -1;

//            //std::map<std::string, int> classNames = database.getClassNames();
//            for (uint i = 0; i < database.getNoFeatures(); ++i)
//            {
//                std::map<std::string, double> classAverages;
//                std::map<std::string, double> classStds;

//                for (auto const &ob : database.getObjects())
//                {
//                    classAverages[ob.getClassName()] += ob.getFeatures()[i];
//                    classStds[ob.getClassName()] += ob.getFeatures()[i] * ob.getFeatures()[i];
//                }

//                std::for_each(database.getClassCounters().begin(), database.getClassCounters().end(), [&](const std::pair<std::string, int> &it)
//                {
//                    classAverages[it.first] /= it.second;
//                    classStds[it.first] = std::sqrt(classStds[it.first] / it.second - classAverages[it.first] * classAverages[it.first]);
//                }
//                );

//                tmp = std::abs(classAverages[ database.getClassNames()[0] ] - classAverages[database.getClassNames()[1]]) / (classStds[database.getClassNames()[0]] + classStds[database.getClassNames()[1]]);

//                if (tmp > FLD)
//                {
//                    FLD = tmp;
//                    max_ind = i;
//                }

//              }

//            ui->FStextBrowserDatabaseInfo->append("max_ind: "  +  QString::number(max_ind) + " " + QString::number(FLD));
//          }
        }
    }
}



void MainWindow::on_FSpushButtonSaveFile_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this,
    tr("Open TextFile"), "D:\\Users\\Krzysiu\\Documents\\Visual Studio 2015\\Projects\\SMPD\\SMPD\\Debug\\", tr("Texts Files (*.txt)"));

        QMessageBox::information(this, "My File", fileName);
        database.save(fileName.toStdString());
}

void MainWindow::on_PpushButtonSelectFolder_clicked()
{
}

void MainWindow::on_CpushButtonOpenFile_clicked()
{

}

void MainWindow::on_CpushButtonSaveFile_clicked()
{

}

void MainWindow::on_CpushButtonTrain_clicked()
{

}

void MainWindow::on_CpushButtonExecute_clicked()
{

}
