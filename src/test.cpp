#include "FinitePowerAlgebras.h"
#include<iostream>

using namespace std;

template<typename GammaType>
auto z(const GammaType& gamma) {
    return gamma(0).conj() * gamma(1) + gamma(3);
}


typedef AltProductRelation<
            CliffordCommRelation<SelfConjRelation<2>>, 
            ExteriorCommRelation<SelfConjRelation<2>>
        > customRelation;

int main(){
    // BaseAlgebra<customRelation> gamma; 
    CliffordAlgebra<3> gamma;
    // auto x = gamma(2) * gamma(1) * gamma(1) + Complex(4., 3.);
    // cout << gamma.anticommutator(gamma(0), gamma(1)) << endl;
    // cout << x << endl;
    // cout << x.norm() << endl;

    CARAlgebra<24> alpha; 
    auto result = alpha.one();
    for (auto j=0; j<12; j++) {
        result = result * alpha(2*j+1) * alpha(2*j);
    }
    cout << result << endl;
    return 0;

}