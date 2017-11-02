#include "vectorcl.h"

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

void printv(std::vector<std::vector<int>> v)
{
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[i].size(); j++) {
            std::cout << " " << v[i][j];
        }
        std::cout << std::endl;
    }
}
