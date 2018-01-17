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

std::vector<int> mv;

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
            mv = maxVector;
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
            mv = maxVector;
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
    ui->CcomboBoxClassifiers->addItem("");
    ui->CcomboBoxClassifiers->addItem("NN");
    ui->CcomboBoxClassifiers->addItem("k-NN");
    ui->CcomboBoxClassifiers->addItem("NM");
    ui->CcomboBoxClassifiers->addItem("k-NM");
    ui->CcomboBoxMethods->addItem("");
    ui->CcomboBoxMethods->addItem("Crossvalidation");
    ui->CcomboBoxMethods->addItem("Bootstrap");
    for (int i = 1; i < 10; i++)
        ui->CcomboBoxK->addItem(QString::number(i));
    ui->CcomboBoxSegs->addItem(QString::number(2));
    ui->CcomboBoxSegs->addItem(QString::number(4));
    ui->CcomboBoxSegs->addItem(QString::number(5));
    ui->CcomboBoxSegs->addItem(QString::number(10));
    train.push_back(0);
    train.push_back(0);
}

void MainWindow::on_CpushButtonSaveFile_clicked()
{

}
int f = 0;
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
    f = 0;
}

std::vector<Object> class1copy;
std::vector<Object> class2copy;

void MainWindow::on_CpushButtonExecute_clicked()
{
    int trainA = train[0];
    int trainB = train[1];
    std::vector<Object> class1;
    std::vector<Object> class2;
    std::vector<Object> all_obj = database.getObjects();
    std::string cur_name = database.getObjects()[0].getClassName();
    if (f == 0){
        for (uint i = 0; i < all_obj.size(); i++) {
            if (all_obj[i].getClassName() == cur_name) {
                class1.push_back(all_obj[i]);
            } else {
                class2.push_back(all_obj[i]);
            }
        }
        random_shuffle(class1.begin(), class1.end());
        random_shuffle(class2.begin(), class2.end());
        class1copy = class1;
        class2copy = class2;
        f = 1;
    } else {
        class1 = class1copy;
        class2 = class2copy;
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
    if (mv.size() > 1) {
        class1ObFeatures = getMatrixFromVector(mv, class1ObFeatures);
        class2ObFeatures = getMatrixFromVector(mv, class2ObFeatures);
        cl1t = getMatrixFromVector(mv, cl1t);
        cl2t = getMatrixFromVector(mv, cl2t);
        for(auto i = mv.begin(); i != mv.end(); ++i) {
            std::cout << *i << " ";
          }
    }
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
        std::vector<std::vector<double>> tA;
        std::vector<std::vector<double>> execA;
        std::vector<std::vector<double>> tB;
        std::vector<std::vector<double>> execB;
        for (int i = 0; i < class1ObFeatures.size(); i++) {
            if (i < trainA) tA.push_back(class1ObFeatures[i]);
            else execA.push_back(class1ObFeatures[i]);
        }
        for (int i = 0; i < class2ObFeatures.size(); i++) {
            if (i < trainB) tB.push_back(class2ObFeatures[i]);
            else execB.push_back(class2ObFeatures[i]);
        }
        double p = NN(tA, execA, tB, execB);
        ui->CtextBrowser->append("NN:\nPercent of correct classifications from "
                                 + QString::number(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB)
                                 + " observations: " + QString::number(p) + "%");
    }
    if (classifier == "NM") {
        std::vector<std::vector<double>> tA;
        std::vector<std::vector<double>> execA;
        std::vector<std::vector<double>> tB;
        std::vector<std::vector<double>> execB;
        for (int i = 0; i < class1ObFeatures.size(); i++) {
            if (i < trainA) tA.push_back(class1ObFeatures[i]);
            else execA.push_back(class1ObFeatures[i]);
        }
        for (int i = 0; i < class2ObFeatures.size(); i++) {
            if (i < trainB) tB.push_back(class2ObFeatures[i]);
            else execB.push_back(class2ObFeatures[i]);
        }
        double p = NM(tA, execA, tB, execB);
        ui->CtextBrowser->append("NM:\nPercent of correct classifications from "
                                 + QString::number(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB)
                                 + " observations: " + QString::number(p) + "%");
    }
    if (classifier == "k-NN") {
        std::vector<std::vector<double>> tA;
        std::vector<std::vector<double>> execA;
        std::vector<std::vector<double>> tB;
        std::vector<std::vector<double>> execB;
        for (int i = 0; i < class1ObFeatures.size(); i++) {
            if (i < trainA) tA.push_back(class1ObFeatures[i]);
            else execA.push_back(class1ObFeatures[i]);
        }
        for (int i = 0; i < class2ObFeatures.size(); i++) {
            if (i < trainB) tB.push_back(class2ObFeatures[i]);
            else execB.push_back(class2ObFeatures[i]);
        }
        double p = kNN(tA, execA, tB, execB, ui->CcomboBoxK->currentText().toInt());
        ui->CtextBrowser->append("k-NN(k=" + QString::number(ui->CcomboBoxK->currentText().toInt()) + "):\nPercent of correct classifications from "
                                 + QString::number(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB)
                                 + " observations: " + QString::number(p) + "%");
    }
    if (classifier == "k-NM") {
        std::vector<std::vector<double>> tA;
        std::vector<std::vector<double>> execA;
        std::vector<std::vector<double>> tB;
        std::vector<std::vector<double>> execB;
        for (int i = 0; i < class1ObFeatures.size(); i++) {
            if (i < trainA) tA.push_back(class1ObFeatures[i]);
            else execA.push_back(class1ObFeatures[i]);
        }
        for (int i = 0; i < class2ObFeatures.size(); i++) {
            if (i < trainB) tB.push_back(class2ObFeatures[i]);
            else execB.push_back(class2ObFeatures[i]);
        }
        double p = kNM(tA, execA, tB, execB, ui->CcomboBoxK->currentText().toInt());
        ui->CtextBrowser->append("k-NM(k=" + QString::number(ui->CcomboBoxK->currentText().toInt()) + "):\nPercent of correct classifications from "
                                 + QString::number(class1ObFeatures.size() + class2ObFeatures.size() - trainA - trainB)
                                 + " observations: " + QString::number(p) + "%");
    }
    if (ui->CcomboBoxMethods->currentText().toStdString() == "Crossvalidation") {
        double p = 0.0;
        for (int i = 0; i < ui->CcomboBoxSegs->currentText().toInt(); i++) {
            int trainAend = (i + 1) * (floor(class1.size()*((100/ui->CcomboBoxSegs->currentText().toInt())/100.0)));
            int trainBend = (i + 1) * (floor(class2.size()*((100/ui->CcomboBoxSegs->currentText().toInt())/100.0)));
            int trainAstart = i * (floor(class1.size()*((100/ui->CcomboBoxSegs->currentText().toInt())/100.0)));
            int trainBstart = i * (floor(class2.size()*((100/ui->CcomboBoxSegs->currentText().toInt())/100.0)));
            if ((i - 1) == ui->CcomboBoxSegs->currentText().toInt()) {
                trainAend = class1ObFeatures.size();
                trainBend = class2ObFeatures.size();
            }
            std::vector<std::vector<double>> tA;
            std::vector<std::vector<double>> execA;
            std::vector<std::vector<double>> tB;
            std::vector<std::vector<double>> execB;
            for (int i = 0; i < class1ObFeatures.size(); i++) {
                if ((i >= trainAstart) && (i < trainAend)) tA.push_back(class1ObFeatures[i]);
                else execA.push_back(class1ObFeatures[i]);
            }
            for (int i = 0; i < class2ObFeatures.size(); i++) {
                if ((i >= trainBstart) && (i < trainBend)) tB.push_back(class2ObFeatures[i]);
                else execB.push_back(class2ObFeatures[i]);
            }
            if (classifier == "NN")
                p = p + NN(tA, execA, tB, execB);
            if (classifier == "k-NN")
                p = p + kNN(tA, execA, tB, execB, ui->CcomboBoxK->currentText().toInt());
            if (classifier == "NM")
                p = p + NM(tA, execA, tB, execB);
            if (classifier == "k-NM")
                p = p + kNM(tA, execA, tB, execB, ui->CcomboBoxK->currentText().toInt());
        }
        ui->CtextBrowser->append("Crossvalidation: " + QString::number(p/ui->CcomboBoxSegs->currentText().toInt()) + "%");
    }
    if (ui->CcomboBoxMethods->currentText().toStdString() == "Bootstrap") {
        std::vector<int> indexesA;
        std::vector<int> indexesB;
        srand(time(0));
        for (int i = 0; i < class1ObFeatures.size(); i++) {
            indexesA.push_back(0 + rand() % (class1ObFeatures.size() - 1));
        }
        for (int i = 0; i < class2ObFeatures.size(); i++) {
            indexesB.push_back(0 + rand() % (class2ObFeatures.size() - 1));
        }
        std::vector<int> tA;
        std::vector<int> tB;
        for (int i = 0; i < indexesA.size(); i++) {
            std::vector<int>::iterator it = std::find(tA.begin(), tA.end(), indexesA[i]);
            if (it == tA.end())
                tA.push_back(indexesA[i]);
        }
        for (int i = 0; i < indexesB.size(); i++) {
            std::vector<int>::iterator it = std::find(tB.begin(), tB.end(), indexesB[i]);
            if (it == tB.end())
                tB.push_back(indexesB[i]);
        }
        std::vector<std::vector<double>> trA;
        std::vector<std::vector<double>> execA;
        std::vector<std::vector<double>> trB;
        std::vector<std::vector<double>> execB;
        for (int i = 0; i < class1ObFeatures.size(); i++) {
            std::vector<int>::iterator it = std::find(tA.begin(), tA.end(), i);
            if (it != tA.end())
                trA.push_back(class1ObFeatures[i]);
            else
                execA.push_back(class1ObFeatures[i]);
        }
        for (int i = 0; i < class2ObFeatures.size(); i++) {
            std::vector<int>::iterator it = std::find(tB.begin(), tB.end(), i);
            if (it != tB.end())
                trB.push_back(class2ObFeatures[i]);
            else
                execB.push_back(class2ObFeatures[i]);
        }
        double p = 0.0;
        if (classifier == "NN")
            p = NN(trA, execA, trB, execB);
        if (classifier == "k-NN")
            p = kNN(trA, execA, trB, execB, ui->CcomboBoxK->currentText().toInt());
        if (classifier == "NM")
            p = NM(trA, execA, trB, execB);
        if (classifier == "k-NM")
            p = kNM(trA, execA, trB, execB, ui->CcomboBoxK->currentText().toInt());
         ui->CtextBrowser->append("Bootstrap: Percent of correct classifications from " +
                                  QString::number(execA.size() + execB.size()) + ": " + QString::number(p) + "%");
    }
}
