#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <QString>

std::vector<std::vector<int>> comb(int N, int K);
float getVectorModule(std::vector<float> v);
std::vector<float> getVectorFromVector(std::vector<int> v, std::vector<float> vv);
std::vector<float> getMatrixMedian(std::vector<std::vector<float>> v);
std::vector<float> getProbabilityVector(int n);
std::vector<float> getVectorDifference(std::vector<float> v1, std::vector<float> v2);
std::vector<std::vector<float>> transponate(std::vector<std::vector<float>> v1);
std::vector<std::vector<float>> setProbability(std::vector<float> prob, std::vector<std::vector<float>> v);
std::vector<std::vector<float>> minusAvg(std::vector<float> v, std::vector<std::vector<float>> m);
std::vector<std::vector<float>> multiplyMatrix(std::vector<std::vector<float>> v1, std::vector<std::vector<float>> v2);
std::vector<std::vector<float>> getMatrixFromVector(std::vector<int> v, std::vector<std::vector<float>> m);
float getMatrixDeterminant(const std::vector<std::vector<float> > &m);
void printv(std::vector<float> v);
void printvv(std::vector<std::vector<float>> v);
std::vector<std::vector<float>> minor(std::vector<std::vector<float>> m, const int &i, const int &j);
QString vectorToString(std::vector<int> v);
