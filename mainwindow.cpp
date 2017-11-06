#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include "vectorcl.h"


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
            std::vector<std::vector<float>> class1ObFeatures;
            std::vector<std::vector<float>> class2ObFeatures;
            for (uint i = 0; i < class1.size(); i++) {
                class1ObFeatures.push_back(class1[i].getFeatures());
            }
            for (uint i = 0; i < class2.size(); i++) {
                class2ObFeatures.push_back(class2[i].getFeatures());
            }
            class1ObFeatures = transponate(class1ObFeatures);
            class2ObFeatures = transponate(class2ObFeatures);
            std::vector<std::vector<float>> cl1MAvg = minusAvg(getMatrixMedian(class1ObFeatures), class1ObFeatures);
            std::vector<std::vector<float>> cl2MAvg = minusAvg(getMatrixMedian(class2ObFeatures), class2ObFeatures);
            std::vector<std::vector<int>> featureCombinations = comb(database.getNoFeatures(), dimension);
            std::map<std::vector<int>, float> ftrs;
            for (int i = 0; i < featureCombinations.size(); i++) {
                std::vector<std::vector<float>> comb1 = getMatrixFromVector(featureCombinations[i], cl1MAvg);
                std::vector<std::vector<float>> comb2 = getMatrixFromVector(featureCombinations[i], cl2MAvg);
                std::vector<std::vector<float>> mm1 = multiplyMatrix(comb1, transponate(comb1), getProbabilityVector(class1.size()));
                std::vector<std::vector<float>> mm2 = multiplyMatrix(comb2, transponate(comb2), getProbabilityVector(class2.size()));
                std::vector<float> median1 = getVectorFromVector(featureCombinations[i], getMatrixMedian(class1ObFeatures));
                std::vector<float> median2 = getVectorFromVector(featureCombinations[i], getMatrixMedian(class2ObFeatures));
                float medianModule = getVectorModule(getVectorDifference(median1, median2));
                float dt1 = getMatrixDeterminant(mm1);
                float dt2 = getMatrixDeterminant(mm2);
                if (dimension == 1) {
                    ftrs[featureCombinations[i]] = medianModule/(sqrt(getMatrixDeterminant(mm1)) + sqrt(getMatrixDeterminant(mm2)));
                } else {
                ftrs[featureCombinations[i]] = medianModule/(getMatrixDeterminant(mm1) + getMatrixDeterminant(mm2));
                }
            }
            float max = 0.0;
            std::pair<std::vector<int>, float> maxp;
            for (auto it = ftrs.begin(); it != ftrs.end(); ++it) {
                if (it->second > max) {
                    max = it->second;
                    maxp = std::make_pair(it->first, it->second);
                }
            }
            QString s = vectorToString(maxp.first);
            ui->FStextBrowserDatabaseInfo->append("max_vector: {"  +  s + "} " + QString::number(maxp.second));
//        {
//            float FLD = 0, tmp;
//            int max_ind = -1;

//            //std::map<std::string, int> classNames = database.getClassNames();
//            for (uint i = 0; i < database.getNoFeatures(); ++i)
//            {
//                std::map<std::string, float> classAverages;
//                std::map<std::string, float> classStds;

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
