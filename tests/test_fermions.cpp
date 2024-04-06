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
}

void derivative_tests() {
    return; 
}

void relation_tests() {
    const uint n = 2;
    DiracAlgebra<n> a; 
    auto l = a.lAlg(); 
    auto r = a.rAlg();
    ComplexDispType d1({
        FieldType(1., 2.), FieldType(5., 6.)
    }), d2({
        FieldType(3., 4.), FieldType(7., 8.)
    }); 
    auto D1 = a.disp(d1), D2=a.disp(d2);
    assert (D1.conj() * a(0) * D1 == a.a(0) + a.disp_vec(d1)[0]); 
    assert (a.disp(add_vec(d1, d2)) == D1 * D2 * a.disp_add(d1, d2)); 
    assert (a.disp(d1).conj() == a.disp(neg_vec(d1))); 
}

int main() {
    majorana_tests();
    fourier_tests();
    derivative_tests();
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