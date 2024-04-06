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

    DiracAlgebra<n> a; 
    auto l = a.lAlg(); 
    auto r = a.rAlg();
    ComplexDispType d1({
        FieldType(1., 2.), FieldType(5., 6.)
    }), d2({
        FieldType(3., 4.), FieldType(7., 8.)
    }); 
    auto D1 = a.disp(d1), D2=a.disp(d2);
    cout << (D1.conj() * a(0) * D1 == a.a(0) + a.disp_vec(d1)[0]) << endl;
    cout << (a.disp(add_vec(d1, d2)) == D1 * D2 * a.disp_add(d1, d2)) << endl;
    cout << (a.disp(d1).conj() == a.disp(neg_vec(d1))) << endl;

    auto c1 = a.coherent(d1), c2=a.coherent(d2); 
    cout << c1 << endl;
}


// auto func = [&zeta, &eta, &alpha](int i) -> auto {
//     return zeta.extL(eta.extR(alpha(i)));
// };