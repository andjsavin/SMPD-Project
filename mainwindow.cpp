#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include "vectorcl.h"
#include "matrixutil.hpp"
#include "invert_matrix.hpp"


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
                            maxVector= vh;
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
    if( ui->FSradioButtonSFS ->isChecked())
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
            double max = 0.0;
            std::vector<int> maxVector;
            std::vector<int> vh;
            std::vector<int> allN;
            for (int i = 0; i < N; i++) {
                allN.push_back(i);
            }
            for (int i = 0; i < K; i++) {
                double tempmax = 0.0;
                maxVector.push_back(-1);
                vh.push_back(-1);
                int j = 0;
                if (i > 0) {
                    for (int x = 0; x < allN.size(); x++) {
                        if (allN[x] == vh[vh.size() -2]) {
                            allN.erase(allN.begin() + x);
                            std::vector<int>(allN).swap(allN);
                        }
                    }
                }
                do {
                    vh[vh.size() - 1] = allN[j];
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
                    if (i == 0) {
                        if (tempmax < medianModule/(sqrt(dt1) + sqrt(dt2)))
                        {
                            tempmax = medianModule/(sqrt(dt1) + sqrt(dt2));
                            maxVector[maxVector.size() - 1] = allN[j];
                        }

                    } else {
                        if (tempmax < medianModule/(dt1 + dt2))
                        {
                            tempmax = medianModule/(dt1 + dt2);
                            maxVector[maxVector.size() - 1] = allN[j];
                        }
                    }
                    j++;
                } while (j < N - i);
                max = tempmax;
                vh[vh.size() - 1] = maxVector[maxVector.size() - 1];
            }
            QString s = vectorToString(maxVector);
            ui->FStextBrowserDatabaseInfo->append("max_vector: {"  +  s + "} " + QString::number(max));
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

std::vector<int> train;

void MainWindow::on_CpushButtonOpenFile_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Open TextFile"), "", tr("Texts Files (*.txt)"));

    if ( !database.load(fileName.toStdString()) )
        QMessageBox::warning(this, "Warning", "File corrupted !!!");
    else
        QMessageBox::information(this, fileName, "File loaded !!!");

    ui->CcomboBoxClassifiers->addItem("k-NN");
    ui->CcomboBoxClassifiers->addItem("k-NM");
    ui->CcomboBoxClassifiers->addItem("NN");
    ui->CcomboBoxClassifiers->addItem("NM");
    for (int i = 2; i < 10; i++)
        ui->CcomboBoxK->addItem(QString::number(i));
    train.push_back(0);
    train.push_back(0);
}

void MainWindow::on_CpushButtonSaveFile_clicked()
{

}

void MainWindow::on_CpushButtonTrain_clicked()
{
    int percent = ui->CplainTextEditTrainingPart->toPlainText().toInt();
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
    train[0] = (floor(class1.size()*(percent/100.0)));
    train[1] = (floor(class2.size()*(percent/100.0)));
}

void MainWindow::on_CpushButtonExecute_clicked()
{
    int trainA = train[0];
    int trainB = train[1];
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
    std::vector<std::vector<double>> cl1t;
    std::vector<std::vector<double>> cl2t;
    for (uint i = 0; i < class1.size(); i++) {
        class1ObFeatures.push_back(class1[i].getFeatures());
        if (i < trainA) {
            cl1t.push_back(class1[i].getFeatures());
        }
    }
    for (uint i = 0; i < class2.size(); i++) {
        class2ObFeatures.push_back(class2[i].getFeatures());
        if (i < trainB) {
            cl2t.push_back(class2[i].getFeatures());
        }
    }
    class1ObFeatures = transponate(class1ObFeatures);
    class2ObFeatures = transponate(class2ObFeatures);
    cl1t = transponate(cl1t);
    cl2t = transponate(cl2t);
    std::vector<double> medianA = getMatrixMedian(cl1t);
    std::vector<double> medianB = getMatrixMedian(cl2t);
    class1ObFeatures = transponate(class1ObFeatures);
    class2ObFeatures = transponate(class2ObFeatures);
    cl1t = transponate(cl1t);
    cl2t = transponate(cl2t);
//    std::vector<std::vector<double>> avgA = minusAvg(medianA, cl1t);
//    std::vector<std::vector<double>> avgB = minusAvg(medianB, cl2t);
//    std::vector<std::vector<double>> mA = multiplyMatrix(avgA, transponate(avgA), getProbabilityVector(trainA));
//    std::vector<std::vector<double>> mB = multiplyMatrix(avgB, transponate(avgB), getProbabilityVector(trainB));
//    bnu::matrix<double> mmA = getMatrix(mA);
//    bnu::matrix<double> mmB = getMatrix(mB);
//    bnu::matrix<double> inverseA(mA.size(), mA.size()), inverseB(mB.size(), mB.size());
//    std::cout << mmA << std::endl;
//    std::cout << mmB << std::endl;
//    InvertMatrix(mmA, inverseA);
//    InvertMatrix(mmB, inverseB);
//    std::cout << inverseA;
//    std::cout << inverseB;
    string classifier = ui->CcomboBoxClassifiers->currentText().toStdString();
    if (classifier == "NN") {
        int correct = 0;
        for (int i = trainA; i < class1ObFeatures.size(); i++) {
            double minA = getDistance(class1ObFeatures[0], class1ObFeatures[i]);
            double minB = getDistance(class2ObFeatures[0], class1ObFeatures[i]);
            for (int j = 1; j < trainA; j++) {
                if (getDistance(class1ObFeatures[j], class1ObFeatures[i]) < minA)
                    minA = getDistance(class1ObFeatures[j], class1ObFeatures[i]);
            }
            for (int j = 1; j < trainB; j++) {
                if (getDistance(class2ObFeatures[j], class1ObFeatures[i]) < minA)
                    minB = getDistance(class2ObFeatures[j], class1ObFeatures[i]);
            }
            if (minA < minB)
                correct++;
        }
        for (int i = trainB; i < class2ObFeatures.size(); i++) {
            double minA = getDistance(class1ObFeatures[0], class2ObFeatures[i]);
            double minB = getDistance(class2ObFeatures[0], class2ObFeatures[i]);
            for (int j = 1; j < trainA; j++) {
                if (getDistance(class1ObFeatures[j], class2ObFeatures[i]) < minA)
                    minA = getDistance(class1ObFeatures[j], class2ObFeatures[i]);
            }
            for (int j = 1; j < trainB; j++) {
                if (getDistance(class2ObFeatures[j], class2ObFeatures[i]) < minA)
                    minB = getDistance(class2ObFeatures[j], class2ObFeatures[i]);
            }
            if (minB < minA)
                correct++;
        }
        double p = (correct*1.0/(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB))*100;
        ui->CtextBrowser->append("NN:\nNumber of correct classifications: "  +  QString::number(correct) + " of "
                                 + QString::number(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB)
                                 + "\n" + QString::number(p) + "%");
    }
    if (classifier == "NM") {
        int correct = 0;
        for (int i = trainA; i < class1ObFeatures.size(); i++) {
            double min = getDistance(medianB, class1ObFeatures[i]);
            if (getDistance(medianA, class1ObFeatures[i]) < min)
                correct++;
        }
        for (int i = trainB; i < class2ObFeatures.size(); i++) {
            double min = getDistance(medianA, class2ObFeatures[i]);
            if (getDistance(medianB, class2ObFeatures[i]) < min)
                correct++;
        }
        double p = (correct*1.0/(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB))*100;
        ui->CtextBrowser->append("NM:\nNumber of correct classifications: "  +  QString::number(correct) + " of "
                                 + QString::number(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB)
                                 + "\n" + QString::number(p) + "%");
    }
    if (classifier == "k-NN") {
        int correct = 0;
        for (int i = trainA; i < class1ObFeatures.size(); i++) {
            std::vector<double> dA;
            std::vector<double> dB;
            for (int j = 0; j < trainA; j++) {
                dA.push_back(getDistance(class1ObFeatures[j], class1ObFeatures[i]));
            }
            for (int j = 0; j < trainB; j++) {
                dB.push_back(getDistance(class2ObFeatures[j], class1ObFeatures[i]));
            }
            sort(dA.begin(), dA.end());
            sort(dB.begin(), dB.end());
            int ai = 0;
            int bi = 0;
            int cb = 0;
            int ca = 0;
            for (int j = 0; j < ui->CcomboBoxK->currentText().toInt(); j++) {
                if (dA[ai] < dB[bi]) {
                    ca++;
                    ai++;
                } else {
                    cb++;
                    bi++;
                }
            }
            if (ca > cb)
                correct++;
        }
        for (int i = trainB; i < class2ObFeatures.size(); i++) {
            std::vector<double> dA;
            std::vector<double> dB;
            for (int j = 0; j < trainA; j++) {
                dA.push_back(getDistance(class1ObFeatures[j], class2ObFeatures[i]));
            }
            for (int j = 0; j < trainB; j++) {
                dB.push_back(getDistance(class2ObFeatures[j], class2ObFeatures[i]));
            }
            sort(dA.begin(), dA.end());
            sort(dB.begin(), dB.end());
            int ai = 0;
            int bi = 0;
            int cb = 0;
            int ca = 0;
            for (int j = 0; j < ui->CcomboBoxK->currentText().toInt(); j++) {
                if (dA[ai] < dB[bi]) {
                    ca++;
                    ai++;
                } else {
                    cb++;
                    bi++;
                }
            }
            if (cb > ca)
                correct++;
        }
        double p = (correct*1.0/(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB))*100;
        ui->CtextBrowser->append("k-NN (k=" + QString::number( ui->CcomboBoxK->currentText().toInt()) +
                                 "):\nNumber of correct classifications: "  +  QString::number(correct) + " of "
                                 + QString::number(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB)
                                 + "\n" + QString::number(p) + "%");
    }
    if (classifier == "k-NM") {
        int correct = 0;
        class1ObFeatures = transponate(class1ObFeatures);
        class2ObFeatures = transponate(class2ObFeatures);
        std::vector<double> mA = getMatrixMedian(class1ObFeatures);
        std::vector<double> mB = getMatrixMedian(class2ObFeatures);
        std::vector<double> mcA = getMatrixMedian(class1ObFeatures);
        std::vector<double> mcB = getMatrixMedian(class2ObFeatures);
        class1ObFeatures = transponate(class1ObFeatures);
        class2ObFeatures = transponate(class2ObFeatures);
        do {
            correct = 0;
            std::vector<std::vector<double>> c1;
            std::vector<std::vector<double>> c2;
            for (uint i = 0; i < class1.size(); i++) {
                if (i < trainA) {
                    c1.push_back(class1[i].getFeatures());
                }
            }
            for (uint i = 0; i < class2.size(); i++) {
                if (i < trainB) {
                    c2.push_back(class2[i].getFeatures());
                }
            }
            for (int i = trainA; i < class1ObFeatures.size(); i++) {
                double min = getDistance(mB, class1ObFeatures[i]);
                if (getDistance(mA, class1ObFeatures[i]) < min) {
                    c1.push_back(class1ObFeatures[i]);
                    correct++;
                } else {
                    c2.push_back(class1ObFeatures[i]);
                }
            }
            for (int i = trainB; i < class2ObFeatures.size(); i++) {
                double min = getDistance(mA, class2ObFeatures[i]);
                if (getDistance(mB, class2ObFeatures[i]) < min) {
                    c2.push_back(class2ObFeatures[i]);
                    correct++;
                } else {
                    c1.push_back(class2ObFeatures[i]);
                }
            }
            c2 = transponate(c2);
            c1 = transponate(c1);
            mcA = getMatrixMedian(c1);
            mcB = getMatrixMedian(c2);
            if (vectCompare(mA, mcA))
                break;
            else {
                mA = mcA;
                mB = mcB;
            }
        } while (true);
        double p = (correct*1.0/(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB))*100;
        ui->CtextBrowser->append("k-NM (k=2):\nNumber of correct classifications: "  +  QString::number(correct) + " of "
                                 + QString::number(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB)
                                 + "\n" + QString::number(p) + "%");
    }
}
