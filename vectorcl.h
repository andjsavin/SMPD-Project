#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <QString>

std::vector<std::vector<int>> comb(int N, int K);
double getVectorModule(std::vector<double> v);
std::vector<double> getVectorFromVector(std::vector<int> v, std::vector<double> vv);
std::vector<double> getMatrixMedian(std::vector<std::vector<double>> v);
std::vector<double> getProbabilityVector(int k);
std::vector<double> getVectorDifference(std::vector<double> v1, std::vector<double> v2);
std::vector<std::vector<double>> transponate(std::vector<std::vector<double>> v1);
std::vector<std::vector<double>> setProbability(std::vector<double> prob, std::vector<std::vector<double>> v);
std::vector<std::vector<double>> minusAvg(std::vector<double> v, std::vector<std::vector<double>> m);
std::vector<std::vector<double>> multiplyMatrix(std::vector<std::vector<double>> v1, std::vector<std::vector<double>> v2, std::vector<double> prob);
std::vector<std::vector<double>> getMatrixFromVector(std::vector<int> v, std::vector<std::vector<double>> m);
double getMatrixDeterminant(const std::vector<std::vector<double> > &m);
void printv(std::vector<double> v);
void printvv(std::vector<std::vector<double>> v);
std::vector<std::vector<double>> minor(std::vector<std::vector<double>> m, const int &i, const int &j);
QString vectorToString(std::vector<int> v);
std::vector<std::vector<double> > mx(double x, std::vector<std::vector<double>> m);
