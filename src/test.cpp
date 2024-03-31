#include "algebras.h"
#include<iostream>

using namespace std;

template<typename GammaType>
auto z(const GammaType& gamma) {
    return gamma(0).conj() * gamma(1) + gamma(3);
}

int main(){
    ExtConjAlg<5> alpha;
    cout << (z(alpha)*z(alpha).conj()) << endl;
    return 0;
}