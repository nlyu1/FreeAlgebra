#include "algebras.h"
#include<iostream>

using namespace std;

int main(){
    CoeffMap coeffs1 = {
        {{0, 1, 0, 0}, Complex(1, 0)},
        {{1, 1, 0, 0}, Complex(2, 0)}
    };
    CoeffMap coeffs2 = {
        {{1, 0, 1, 0}, Complex(3, 0)},
        {{0, 0, 1, 1}, Complex(4, 0)}
    };
    auto x = ExtElm(coeffs1, 4);
    auto y = ExtElm(coeffs2, 4);
    // cout << prettyPrint(x) << endl 
    // << prettyPrint(y) << endl 
    // << prettyPrint(x + (y* Complex(-1, 0))) << endl;
    // cout << prettyPrint(x - y) << endl;

    // cout << prettyPrint(x * y) << endl;

    // CoeffMap c3 = {
    //     {{0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1}, Complex(1., 0.)},
    // };
    // CoeffMap c4 = {
    //     {{1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0}, Complex(1., 0.)},
    // };
    // auto x1 = DiracElm(c3, 32), y1=DiracElm(c4, 32);
    // cout << prettyPrint(x1 * y1) << endl;

    // CoeffMap c5 = {{{1, 1, 1, 1, 1, 0}, Complex(2., 3.)}};
    // auto x2 = DiracElm(c5, 6);
    // // cout << prettyPrint(x2) << endl;
    // cout << prettyPrint(x2) << endl;
    // cout << prettyPrint(x2.conj()) << endl;

    ExtConjAlg alpha(4);
    auto j = alpha.a(0)*alpha.a(2) + alpha.a(1)*alpha.a(3) + alpha.a(3).conj()* alpha.a(2).conj(); 
    cout << ((alpha.log(alpha.exp(j)) - j) == 0) << endl;
    return 0;
}