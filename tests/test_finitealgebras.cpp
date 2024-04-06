#include "utils.h"
#include "BaseAlgebraRelations.h"
#include "FinitePowerAlgebras.h"
#include<iostream>
#include<cassert>
using namespace std; 

void pow_tests() {
    ExtFreeConjAlgebra<2> eta;
    CliffordAlgebra<2> gamma;
    CARAlgebra<2> delta;
    assert (eta(0).pow(0) == 1.);
    assert (eta(0).pow(2) == 0.);
    assert (gamma(0).pow(0) == 1.);
    assert (gamma(0).pow(2) == 1.);
    assert (delta(0).pow(0) == 1.);
    assert (delta(0).pow(2) == 0.);
}

void anticomm_tests() {
    ExtFreeConjAlgebra<2> eta;
    CliffordAlgebra<2> gamma;
    CARAlgebra<2> delta;
    assert (eta.anticommutator(eta(0), eta(1)) == 0.);
    assert (delta.anticommutator(delta(0), delta(1)) == 1.);
    assert (gamma.anticommutator(gamma(0), gamma(1)) == 0.);
}

void exterior_tests() {
    const uint n = 4;
    ExtFreeConjAlgebra<n> b; 
    auto y = b(0) * b(0).conj() * b(2) * b(2).conj(); 
    assert (b.td(y) == 1);
}

int main() {
    pow_tests();
    anticomm_tests();
    exterior_tests();
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