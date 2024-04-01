#include "FinitePowerAlgebras.h"
#include<iostream>

using namespace std;

template<typename GammaType>
auto z(const GammaType& gamma) {
    return gamma(0).conj() * gamma(1) + gamma(3);
}

typedef CliffordAlgebra<3> LAlg;
typedef CARAlgebra<4> RAlg; 

int main(){
    LAlg alpha;
    RAlg beta;
    ProductAlgebra<LAlg, RAlg, decltype(PROD_ANTICOMMUTE)> gamma; 
    // cout << alpha(0) << endl;
    // cout << beta(0) << endl;
    // cout << gamma.extR(alpha(0)) << gamma.extL(beta(0)) << endl;
    // cout << gamma.anticommutator(gamma.extR(alpha(0)), gamma.extL(beta(0))) << endl; 
    // cout << (gamma.projL(gamma.extR(alpha(0)) + gamma.extL(beta(0)))) << endl;
    cout << gamma.extL(beta(0)) << endl;
    cout << (gamma.projL(gamma.extL(beta(0)))) << endl;
    cout << (gamma.projR(gamma.extR(alpha(0))) == beta.zero()) << endl;
    return 0;

}