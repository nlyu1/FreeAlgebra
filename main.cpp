#include "utils.h"
#include "BaseAlgebraRelations.h"
#include "FinitePowerAlgebras.h"
#include "Automorphism.h"
#include "ProductPowerAlgebra.h"
#include "Fermions.h"
#include<iostream>
#include <cassert>
#include<cmath>
using namespace std; 

int main() {
    const uint n = 2;

    DiracAlgebra<n> a; 
    // X = exp(i pi / 2 (1 - (a + adag)))
    auto coeffs = RealCoeffMap({
        {{0, 1}, 100.}, 
        {{0, 1, 2, 3}, -1.},
        {{1, 2, 3}, 2.}
    }); 
    
    auto rho = a.pure(coeffs);
    cout << rho << endl;
    // 
    cout << a.tr(rho * rho.log()) << endl;
    // cout << "Physical trace: " << a.tr(rho) << endl;
}

// Todo: multiplication caching, characteristic function