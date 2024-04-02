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

int main() {
    pow_tests();
    anticomm_tests();
}