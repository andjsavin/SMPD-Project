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

double getVectorModule(std::vector<double> v)
{
    double m =0.0;
    for (int i = 0; i < v.size(); i++) {
        m += v[i]*v[i];
    }
    return sqrt(m);
}

std::vector<double> getMatrixMedian(std::vector<std::vector<double> > v)
{
    std::vector<double> mm;
    for (int i = 0; i < v.size(); i++) {
        double h = 0.0;
        for (int j = 0; j < v[i].size(); j++) {
            h += v[i][j];
        }
        mm.push_back(h/v[i].size());
    }
    return mm;
}

std::vector<double> getVectorDifference(std::vector<double> v1, std::vector<double> v2)
{
    std::vector<double> v;
    for (int i = 0; i < v1.size(); i++) {
        double h = v1[i] - v2[i];
        v.push_back(h);
    }
    return v;
}

std::vector<std::vector<double>> getMatrixFromVector(std::vector<int> v, std::vector<std::vector<double>> m)
{
    std::vector<std::vector<double>> newm;
    for (int i = 0; i < m.size(); i++) {
        if (std::find(v.begin(), v.end(), i) != v.end()) {
            std::vector<double> h;
            for (int j = 0; j < m[i].size(); j++) {
                h.push_back(m[i][j]);
            }
            newm.push_back(h);
        }
    }
    return newm;
}

std::vector<std::vector<double>> minor(std::vector<std::vector<double>> m, const int &i, const int &j) {
    m.erase(m.begin() + i);
    for (auto &a_m : m) {
        a_m.erase(a_m.begin() + j);
    }
    return m;
}

double getMatrixDeterminant(const std::vector<std::vector<double>> &m)
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
    double summ = 0.0;
    int j = 0;
    for (auto &m_0j: m[0]) {
        summ += m_0j*signum*getMatrixDeterminant(minor(m, 0, j));
        signum *= -1;
        j++;
    }
    return summ;
}

std::vector<std::vector<double>> transponate(std::vector<std::vector<double>> v1)
{
    double temp;
    std::vector<std::vector<double>> vh;
    for (int i = 0; i < v1[0].size(); i++) {
        std::vector<double> vf;
        for (int j = 0; j < v1.size(); j++) {
            vf.push_back(v1[j][i]);
        }
        vh.push_back(vf);
    }
    return vh;
}

void printv(std::vector<double> v) {
    for (int i = 0; i < v.size(); i++) {
        std::cout << " " << v[i];
    }
    std::cout << std::endl;
}

void printvv(std::vector<std::vector<double>> v) {
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[i].size(); j++) {
            std::cout << " " << v[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

std::vector<double> getProbabilityVector(int k)
{
    std::vector<double> v;
    for (int i = 0; i < k; i++) {
        v.push_back(1.0/k);
    }
    return v;
}

std::vector<std::vector<double>> minusAvg(std::vector<double> v, std::vector<std::vector<double>> m)
{
    std::vector<std::vector<double>> vh;
    for (int i = 0; i < m.size(); i++) {
        std::vector<double> vv;
        for (int j = 0; j < m[i].size(); j++) {
            vv.push_back(m[i][j] - v[i]);
        }
        vh.push_back(vv);
    }
    return vh;
}

std::vector<std::vector<double>> multiplyMatrix(std::vector<std::vector<double>> v1, std::vector<std::vector<double>> v2, std::vector<double> prob)
{
    std::vector<std::vector<double>> vh;
    for (int i = 0; i < v1.size(); i++) {
        std::vector<double> vv;
        for (int j = 0; j < v2[0].size(); j++) {
            double temp = 0;
            for (int inner = 0; inner < v1[i].size(); inner++) {
                temp += v1[i][inner]*v2[inner][j] * prob[inner];
            }
            vv.push_back(temp);
        }
        vh.push_back(vv);
    }
    return vh;
}

std::vector<double> getVectorFromVector(std::vector<int> v, std::vector<double> vv)
{
    std::vector<double> vh;
    for (int i = 0; i < vv.size(); i++) {
        if (std::find(v.begin(), v.end(), i) != v.end())
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

std::vector<std::vector<double>> mx(double x, std::vector<std::vector<double>> m)
{
    std::vector<std::vector<double>> vh;
    for (int i = 0; i < m.size(); i++){
        std::vector<double> vv;
        for (int j = 0; j < m[i].size(); j++) {
            vv.push_back(m[i][j]/x);
        }
        vh.push_back(vv);
    }
    return vh;
}

double getDistance(std::vector<double> v1, std::vector<double> v2) {
    double res = 0.0;
    for (int i = 0; i < v1.size(); i++) {
        res += (v1[i] - v2[i])*(v1[i] - v2[i]);
    }
    res = sqrt(res);
    return res;
}

bool vectCompare(std::vector<double> v1, std::vector<double> v2) {
    for (int i = 0; i < v1.size(); i++) {
        if (v1[i] != v2[i])
            return false;
    }
    return true;
}

double NN(std::vector<std::vector<double>> trainA, std::vector<std::vector<double>> execA,
          std::vector<std::vector<double>> trainB, std::vector<std::vector<double>> execB) {
    int correct = 0;
    for (int i = 0; i < execA.size(); i++) {
        double minA = getDistance(trainA[0], execA[i]);
        double minB = getDistance(trainB[0], execA[i]);
        for (int j = 1; j < trainA.size(); j++) {
            if (getDistance(trainA[j], execA[i]) < minA)
                minA = getDistance(trainA[j], execA[i]);
        }
        for (int j = 1; j < trainB.size(); j++) {
            if (getDistance(trainB[j], execA[i]) < minB)
                minB = getDistance(trainB[j], execA[i]);
        }
        if (minA < minB)
            correct++;
    }
    for (int i = 0; i < execB.size(); i++) {
        double minA = getDistance(trainA[0], execB[i]);
        double minB = getDistance(trainB[0], execB[i]);
        for (int j = 1; j < trainA.size(); j++) {
            if (getDistance(trainA[j], execB[i]) < minA)
                minA = getDistance(trainA[j], execB[i]);
        }
        for (int j = 1; j < trainB.size(); j++) {
            if (getDistance(trainB[j], execB[i]) < minB)
                minB = getDistance(trainB[j], execB[i]);
        }
        if (minB < minA)
            correct++;
    }
    return (correct*1.0/(execA.size() + execB.size()))*100;
}

double NM(std::vector<std::vector<double>> trainA, std::vector<std::vector<double>> execA,
          std::vector<std::vector<double>> trainB, std::vector<std::vector<double>> execB) {
    int correct = 0;
    std::vector<double> medianA = getMatrixMedian(transponate(trainA));
    std::vector<double> medianB = getMatrixMedian(transponate(trainB));
    for (int i = 0; i < execA.size(); i++) {
        double min = getDistance(medianB, execA[i]);
        if (getDistance(medianA, execA[i]) < min)
            correct++;
    }
    for (int i = 0; i < execB.size(); i++) {
        double min = getDistance(medianA, execB[i]);
        if (getDistance(medianB, execB[i]) < min)
            correct++;
    }
    return (correct*1.0/(execA.size() + execB.size()))*100;
}

bool matrix_comp(std::vector<double> m1, std::vector<double> m2)
{
    if (m1.size() != m2.size()) return false;
    for (int i = 0; i < m1.size(); i++) {
        if (m1[i] != m2[i]) return false;
    }
    return true;
}

double kNN(std::vector<std::vector<double>> trainA, std::vector<std::vector<double>> execA,
          std::vector<std::vector<double>> trainB, std::vector<std::vector<double>> execB, int k) {
    int correct = 0;
    for (int i = 0; i < execA.size(); i++) {
        std::vector<double> dA;
        std::vector<double> dB;
        for (int j = 0; j < trainA.size(); j++) {
            dA.push_back(getDistance(trainA[j], execA[i]));
        }
        for (int j = 0; j < trainB.size(); j++) {
            dB.push_back(getDistance(trainB[j], execA[i]));
        }
        sort(dA.begin(), dA.end());
        sort(dB.begin(), dB.end());
        int ai = 0;
        int bi = 0;
        int cb = 0;
        int ca = 0;
        for (int j = 0; j < k; j++) {
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
    for (int i = 0; i < execB.size(); i++) {
        std::vector<double> dA;
        std::vector<double> dB;
        for (int j = 0; j < trainA.size(); j++) {
            dA.push_back(getDistance(trainA[j], execB[i]));
        }
        for (int j = 0; j < trainB.size(); j++) {
            dB.push_back(getDistance(trainB[j], execB[i]));
        }
        sort(dA.begin(), dA.end());
        sort(dB.begin(), dB.end());
        int ai = 0;
        int bi = 0;
        int cb = 0;
        int ca = 0;
        for (int j = 0; j < k; j++) {
            if (dB[bi] < dA[ai]) {
                cb++;
                bi++;
            } else {
                ca++;
                ai++;
            }
        }
        if (cb > ca)
            correct++;
    }
    return (correct*1.0/(execA.size() + execB.size()))*100;
}

double kNM(std::vector<std::vector<double>> trainA, std::vector<std::vector<double>> execA,
           std::vector<std::vector<double>> trainB, std::vector<std::vector<double>> execB, int k) {
    std::vector<std::vector<double>> meansA;
    std::vector<std::vector<double>> meansB;
    for (int i = 0; i < k; i++) {
        meansA.push_back(trainA[i]);
        meansB.push_back(trainB[i]);
    }
    do {
        std::map<int, std::vector<std::vector<double>>> elementsA;
        std::map<int, std::vector<std::vector<double>>> elementsB;
        for (int i = 0; i < k; i++) {
            elementsA.insert(std::pair<int, std::vector<std::vector<double>>>(i, std::vector<std::vector<double>>()));
            elementsB.insert(std::pair<int, std::vector<std::vector<double>>>(i, std::vector<std::vector<double>>()));
        }
        for (int i = 0; i < trainA.size(); i++) {
            double min = getDistance(trainA[i], meansA[0]);
            int id = 0;
            for (int j = 1; j < meansA.size(); j++) {
                if (getDistance(trainA[i], meansA[j]) < min) {
                    id = j;
                    min = getDistance(trainA[i], meansA[j]);
                }
            }
            std::map<int, std::vector<std::vector<double>>>::iterator it = elementsA.find(id);
            if (it != elementsA.end())
                it->second.push_back(trainA[i]);
        }
        for (int i = 0; i < trainB.size(); i++) {
            double min = getDistance(trainB[i], meansB[0]);
            int id = 0;
            for (int j = 1; j < meansB.size(); j++) {
                if (getDistance(trainB[i], meansB[j]) < min) {
                    id = j;
                    min = getDistance(trainB[i], meansB[j]);
                }
            }
            std::map<int, std::vector<std::vector<double>>>::iterator it = elementsB.find(id);
            if (it != elementsB.end())
                it->second.push_back(trainB[i]);
        }
        std::vector<std::vector<double>> new_meansA;
        std::vector<std::vector<double>> new_meansB;
        for (int i = 0; i < meansA.size(); i++) {
            std::map<int, std::vector<std::vector<double>>>::iterator itA = elementsA.find(i);
            std::map<int, std::vector<std::vector<double>>>::iterator itB = elementsB.find(i);
            itA->second = transponate(itA->second);
            itB->second = transponate(itB->second);
            new_meansA.push_back(getMatrixMedian(itA->second));
            new_meansB.push_back(getMatrixMedian(itB->second));
        }
        int flag = 0;
        for (int i = 0; i < meansA.size(); i++) {
            if (matrix_comp(new_meansA[i], meansA[i]) == false)
            {
                flag = 1;
                meansA[i] = new_meansA[i];
            }
            if (matrix_comp(new_meansB[i], meansB[i]) == false) {
                flag = 1;
                meansB[i] = new_meansB[i];
            }
        }
        if (flag == 0) break;
    } while (true);
    int correct = 0;
    for (int i = 0; i < execA.size(); i++) {
        std::vector<double> dA;
        std::vector<double> dB;
        for (int j = 0; j < meansA.size(); j++) {
            dA.push_back(getDistance(meansA[j], execA[i]));
        }
        for (int j = 0; j < meansB.size(); j++) {
            dB.push_back(getDistance(meansB[j], execA[i]));
        }
        sort(dA.begin(), dA.end());
        sort(dB.begin(), dB.end());
        int ai = 0;
        int bi = 0;
        int cb = 0;
        int ca = 0;
        for (int j = 0; j < k; j++) {
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
    for (int i = 0; i < execB.size(); i++) {
        std::vector<double> dA;
        std::vector<double> dB;
        for (int j = 0; j < meansA.size(); j++) {
            dA.push_back(getDistance(meansA[j], execB[i]));
        }
        for (int j = 0; j < meansB.size(); j++) {
            dB.push_back(getDistance(meansB[j], execB[i]));
        }
        sort(dA.begin(), dA.end());
        sort(dB.begin(), dB.end());
        int ai = 0;
        int bi = 0;
        int cb = 0;
        int ca = 0;
        for (int j = 0; j < k; j++) {
            if (dB[bi] < dA[ai]) {
                cb++;
                bi++;
            } else {
                ca++;
                ai++;
            }
        }
        if (cb > ca)
            correct++;
    }
    return (correct*1.0/(execA.size() + execB.size()))*100;
}
