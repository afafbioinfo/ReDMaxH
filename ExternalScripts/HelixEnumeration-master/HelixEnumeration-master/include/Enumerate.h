#ifndef ENUMERATE_H
#define ENUMERATE_H

#include <vector>
#include <string>

namespace Enumerate {

struct HelixClass {
    int i;
    int j;
    int k;

    HelixClass(int i, int j, int k) 
        : i(i), j(j), k(k) {}
};

std::vector<HelixClass> EnumerateHelices(
        const std::string& sequence,
        const int min_k = 3,
        const int hairpin_length = 3);

} //Enumerate namespace

#endif
