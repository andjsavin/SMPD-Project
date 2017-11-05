#include "vectorcl.h"
#include <math.h>
#include <algorithm>

long double fact(int N) //factorial
{
    if (N < 0) {
        return 0;
    }
    if (N == 0) {
        return 1;
    } else {
        return N * fact(N - 1);
    }
}

std::vector<std::vector<int> > comb(int N, int K)
{
    std::string bitmask(K, 1); // We make a bitmask (string) which has K (number of features) first 1
    bitmask.resize(N, 0); // Then we resize it to fill leftover places with 0s (N-K)
    std::vector<std::vector<int>> vc; // Vector of vectors
    //store integers in vector and permute bitmask
//    int c = 1; //for testing
    do {
        std::vector<int> vh; //helping vector
//        std::cout << c << " | ";
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i]) { // if we have 1 in the N bit of bitmask we store it in subvector
//                std::cout << " " << i; //for testing purposes
                vh.push_back(i);
            }
        }
//        c++;//for testing
//        std::cout << std::endl; //for testing purposes
        vc.push_back(vh); // Store subvector in vector of vectors
    } while (std::prev_permutation(bitmask.begin(), bitmask.end())); // lexicographilly permute bitmask
    //for example if we had [1 1 0 0] we will have [1 0 1 0]
    return vc;
}

float getVectorModule(std::vector<float> v)
{
    float m =0.0;
    for (int i = 0; i < v.size(); i++) {
        m += v[i]*v[i];
    }
    return sqrt(m);
}

std::vector<float> getMatrixMedian(std::vector<std::vector<float> > v)
{
    std::vector<float> mm;
    for (int i = 0; i < v.size(); i++) {
        float h = 0.0;
        for (int j = 0; j < v[i].size(); j++) {
            h += v[i][j];
        }
        mm.push_back(h/v[i].size());
    }
    return mm;
}

std::vector<float> getVectorDifference(std::vector<float> v1, std::vector<float> v2)
{
    std::vector<float> v;
    for (int i = 0; i < v1.size(); i++) {
        float h = v1[i] - v2[i];
        v.push_back(h);
    }
    return v;
}

std::vector<std::vector<float>> getMatrixFromVector(std::vector<int> v, std::vector<std::vector<float>> m)
{
    std::vector<std::vector<float>> newm;
    for (int i = 0; i < m.size(); i++) {
        if (std::binary_search(v.begin(), v.end(), i)) {
            std::vector<float> h;
            for (int j = 0; j < m[i].size(); j++) {
                h.push_back(m[i][j]);
            }
            newm.push_back(h);
        }
    }
    return newm;
}

std::vector<std::vector<float>> minor(std::vector<std::vector<float>> m, const int &i, const int &j) {
    m.erase(m.begin() + i);
    for (auto &a_m : m) {
        a_m.erase(a_m.begin() + j);
    }
    return m;
}

float getMatrixDeterminant(const std::vector<std::vector<float>> &m)
{
    int k = m.size();
    int l = m[0].size();
    if (k != l) {
        return NAN;
    }
    if (k == 1) {
        return m[0][0];
    }
    int signum = 1;
    float summ = 0.0;
    int j = 0;
    for (auto &m_0j: m[0]) {
        summ += m_0j*signum*getMatrixDeterminant(minor(m, 0, j));
        signum *= -1;
        j++;
    }
    return summ;
}

std::vector<std::vector<float>> transponate(std::vector<std::vector<float>> v1)
{
    float temp;
    std::vector<std::vector<float>> vh;
    for (int i = 0; i < v1[0].size(); i++) {
        std::vector<float> vf;
        for (int j = 0; j < v1.size(); j++) {
            vf.push_back(v1[j][i]);
        }
        vh.push_back(vf);
    }
    return vh;
}

void printv(std::vector<float> v) {
    for (int i = 0; i < v.size(); i++) {
        std::cout << " " << v[i];
    }
    std::cout << std::endl;
}

void printvv(std::vector<std::vector<float>> v) {
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[i].size(); j++) {
            std::cout << " " << v[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

std::vector<float> getProbabilityVector(int k)
{
    std::vector<float> v;
    for (int i = 0; i < k; i++) {
        v.push_back(1.0/k);
    }
    return v;
}

std::vector<std::vector<float>> minusAvg(std::vector<float> v, std::vector<std::vector<float>> m)
{
    std::vector<std::vector<float>> vh;
    for (int i = 0; i < m.size(); i++) {
        std::vector<float> vv;
        for (int j = 0; j < m[i].size(); j++) {
            vv.push_back(m[i][j] - v[i]);
        }
        vh.push_back(vv);
    }
    return vh;
}

std::vector<std::vector<float>> multiplyMatrix(std::vector<std::vector<float>> v1, std::vector<std::vector<float>> v2, std::vector<float> prob)
{
    std::vector<std::vector<float>> vh;
    for (int i = 0; i < v1.size(); i++) {
        std::vector<float> vv;
        for (int j = 0; j < v2[0].size(); j++) {
            float temp = 0;
            for (int inner = 0; inner < v1[i].size(); inner++) {
                temp += v1[i][inner]*v2[inner][j] * prob[inner];
            }
            vv.push_back(temp);
        }
        vh.push_back(vv);
    }
    return vh;
}

std::vector<float> getVectorFromVector(std::vector<int> v, std::vector<float> vv)
{
    std::vector<float> vh;
    for (int i = 0; i < vv.size(); i++) {
        if (std::binary_search(v.begin(), v.end(), i))
            vh.push_back(vv[i]);
    }
    return vh;
}

QString vectorToString(std::vector<int> v)
{
    QString s = " ";
    for (int i = 0; i < v.size(); i++) {
        s += QString::number(v[i]) + " ";
    }
    return s;
}

std::vector<std::vector<float>> mx(float x, std::vector<std::vector<float>> m)
{
    std::vector<std::vector<float>> vh;
    for (int i = 0; i < m.size(); i++){
        std::vector<float> vv;
        for (int j = 0; j < m[i].size(); j++) {
            vv.push_back(m[i][j]/x);
        }
        vh.push_back(vv);
    }
    return vh;
}
