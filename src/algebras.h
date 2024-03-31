#ifndef ALGEBRAS_H
#define ALGEBRAS_H
#include "element.h"
using namespace std; 

// Algebraic relations for exterior algebra elements with conjugation
// Canonical ordering (a1, a1*, a2, a2*, ...)
struct ExtConjRelation {
    CoeffMap commute(uint i, uint j) const {
        CoeffMap result;
        if (i < j) {
            throw(std::invalid_argument("Expects non-canonical ordering"));
        }
        if (i > j) {
            result[{j, i}] = Complex(-1., 0.);
        }
        return result; 
    }
    // Specifies which element to conjugate to
    uint conj(uint i) const {
        return (i % 2 == 0) ? i+1 : i-1;
    }
};

// Dirac commutation relation (CAR)
// Canonical ordering (a1, a1*, a2, a2*, ...)
struct DirRelation {
    CoeffMap commute(uint i, uint j) const {
        CoeffMap result;
        if (i < j) {
            throw(std::invalid_argument("Expects non-canonical ordering"));
        }
        if (i > j) {
            if (i % 2 == 1 && j == i - 1) {
                result[{}] = Complex(1., 0.);
            }   
            result[{j, i}] = Complex(-1., 0.);
        }
        return result; 
    }

    uint conj(uint i) const {
        return (i % 2 == 0) ? i+1 : i-1;
    }
};

typedef AlgebraElement<DirRelation> DiracElm;
typedef AlgebraElement<ExtConjRelation> ExtElm;

// Exterior conjugation algebra
class ExtConjAlg {
public:
    uint n; // Number of modes
    
    ExtConjAlg(uint n): n(n) {};

    ExtElm a(uint i) {
        assert(i >= 0 && i<n); 
        KeyType coeffs({2*i});
        return ExtElm({{gen_to_power_repr(coeffs, 2*n), 1}}, 2*n);
    }

    ExtElm exp(ExtElm a) {

    }
};
#endif