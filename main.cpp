#include "FinitePowerAlgebras.h"
#include "Automorphism.h"
#include "Fermions.h"
#include <torch/torch.h>
#include<iostream>

using namespace std;

template<typename GammaType>
auto z(const GammaType& gamma) {
    return gamma(0).conj() * gamma(1) + gamma(3);
}

typedef CliffordAlgebra<3> LAlg;
typedef CARAlgebra<4> RAlg; 

int main(){
    // LAlg alpha;
    // RAlg beta;
    // ProductAlgebra<LAlg, RAlg, decltype(PROD_ANTICOMMUTE)> gamma; 
    // cout << alpha(0) << endl;
    // cout << beta(0) << endl;
    // cout << gamma.extR(alpha(0)) << gamma.extL(beta(0)) << endl;
    // cout << gamma.anticommutator(gamma.extR(alpha(0)), gamma.extL(beta(0))) << endl; 
    // cout << (gamma.projL(gamma.extR(alpha(0)) + gamma.extL(beta(0)))) << endl;
    // cout << gamma.extL(beta(0)) << endl;
    // cout << (gamma.projL(gamma.extL(beta(0)))) << endl;
    // cout << (gamma.projR(gamma.extR(alpha(0))) == beta.zero()) << endl;

    const uint n=3;
    MajoranaAlgebra<n> gamma;
    // ProductAlgebra<MajoranaAlgebra<n>, MajoranaAlgebra<n>, decltype(PROD_ANTICOMMUTE)> gamma_sq; 

    RealCoeffMap coeffs = {
        {{}, 3.},
        {{0, 1}, 1.}, 
        {{1, 2}, 2.},
        {{2, 3}, 3.}
    };

    // Gaussian hamiltonian 
    auto rho = gamma.density(coeffs);
    auto U = gamma.unitary(coeffs);
    // cout << "New   " << endl;
    // cout << gamma.moments(rho) << endl;
    // cout << rho << endl;
    // cout << gamma_sq.kron(rho, rho) << endl;
    // cout << gamma_sq.kron(rho, rho).tr() << endl;

    // const uint n=3; 
    // CoeffMap coeffs = {
    //     {{}, Complex(3.)},
    //     {{0, 1}, Complex(1.)}, 
    //     {{1, 2}, Complex(2.)},
    
    //     {{2, 3}, Complex(3.)}
    // };
    // MajoranaAlgebra<n> gamma; 
    // cout << gamma(0).conj() << endl;

    DiracMajoranaBijection<6> R;
    CARAlgebra<6> a;
    ImageRelation<CanonicalAnticommRelation<6>, DiracMajoranaBijection<6>> r;
}