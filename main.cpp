#include "utils.h"
#include "BaseAlgebraRelations.h"
#include "FinitePowerAlgebras.h"
#include "Automorphism.h"
#include "ProductPowerAlgebra.h"
#include "Fermions.h"
#include<iostream>
#include <cassert>
using namespace std; 

int main() {
    const uint n = 2;
    MajoranaAlgebra<n> c;
    // auto rho_coeffs = RealCoeffMap({
    //     {{0, 2, 4, 6}, 1},
    //     {{0, 1, 2, 3}, 1},
    //     {{4, 5, 6, 7}, 1},
    // });
    auto rho_coeffs = RealCoeffMap({
        // {{0, 3}, 1},
        // {{1, 2}, 1},
        // {{2, 3}, 1},
    });
    auto rho = c.pure(rho_coeffs);
    auto xi = c.F(rho); 
    cout << xi << endl;
    auto sigma = c.iF(xi);
    cout << sigma << endl;

    // auto xi = c.F(rho);
    // auto xi_ = c.F_(rho);
    // cout << "Rho " << rho << endl;
    // cout << "xi " << xi << endl;
    // cout << xi_ << endl;
    // cout << xi - xi_ << endl;
}

// Test algebraic relations 
// int main() {
//     const uint n = 2;
//     MajoranaAlgebra<n> c;
//     // auto rho_coeffs = RealCoeffMap({
//     //     {{0, 2, 4, 6}, 1},
//     //     {{0, 1, 2, 3}, 1},
//     //     {{4, 5, 6, 7}, 1},
//     // });
//     auto rho_coeffs = RealCoeffMap({
//         {{0, 3}, 1},
//         {{1, 2}, 1},
//         {{2, 3}, 1},
//     });
//     auto rho = c.pure(rho_coeffs);
//     auto xi = c.F(rho);
//     cout << xi << endl;

//     auto zeta = xi.log();
//     cout << "Exponential consistency " << (zeta.exp() - xi).norm() << endl;
//     cout << zeta << zeta.isreal() << endl;
// }


// If commute: then moment (additional minus sign)
// If non-commute, then no additional minus sign: only scaling