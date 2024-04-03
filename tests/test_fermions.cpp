#include "utils.h"
#include "BaseAlgebraRelations.h"
#include "FinitePowerAlgebras.h"
#include "Automorphism.h"
#include "ProductPowerAlgebra.h"
#include "Fermions.h"
#include<iostream>
#include<cassert>
using namespace std; 

void majorana_tests() {
    MajoranaAlgebra<4> c;
    auto rho_zero = c.zero_state();
    assert ((rho_zero.pow(2) - rho_zero).norm() == 0);
    assert (rho_zero.tr() == 1.);
}

void fourier_tests() {
    // Assert the behavior of fourier transform: preserves coefficients
    const uint n = 3;
    MajoranaAlgebra<n> c;
    auto rho_coeffs = RealCoeffMap({
        {{0, 1}, 1},
        {{2, 3}, 1},
        {{4, 5}, 1},
    });
    auto rho = c.pure_gaussian(rho_coeffs);
    auto xi = c.F(rho);
    for (auto p: xi.coeffs) {
        assert(p.second == rho.coeffs.at(p.first));
    }
    // Assert the behavior of 
}

int main() {
    majorana_tests();
}



// int main() {
//     // Create a 2x2 complex tensor
//     auto real = torch::tensor({{1., 2.}, 
//                                 {3., 4.}});
//     auto imag = torch::tensor({{5., 6.}, 
//                                 {7., 8.}});
//     // real.set_requires_grad(true);
//     // imag.set_requires_grad(true);
//     auto omega = torch::cat({real.unsqueeze(-1), imag.unsqueeze(-1)}, -1);
//     omega = torch::view_as_complex(omega);
//     omega.set_requires_grad(true);
//     cout << omega.sizes() << " " << omega.scalar_type() << endl;
//     auto beta = FieldType(omega[0][0]);

//     beta.tensor().abs().backward();
//     cout << omega.grad() << endl;
//     return 0;
// }