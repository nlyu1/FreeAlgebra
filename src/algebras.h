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
    uint n; // Number of modes, there are 2*n generators in total
    
    ExtConjAlg(uint n): n(n) {};

    ExtElm a(uint i) {
        assert(i >= 0 && i<n); 
        KeyType coeffs({2*i});
        return ExtElm({{gen_to_power_repr(coeffs, 2*n), 1}}, 2*n);
    }

    // Grassmann exponentiation: only need to compute a few terms, yeah
    ExtElm exp(const ExtElm& a) {
        auto power = a.one(), result=a.one(); // Start with the multiplicative identity;
        uint niters = min(static_cast<uint>(a.coeffs.size()), 2*n);
        for (uint j=1; j<=niters; j++) {
            power = power * a;
            result.add_(power * Complex(1./static_cast<double>(factorial(j)), 0));
        }
        return result;
    }

    // Taylor expansion definition of logarithm. 
    ExtElm log(const ExtElm& a) {
        auto b = a - Complex(1.);
        auto power = a.one(), result = a.zero();
        uint niters = min(static_cast<uint>(b.coeffs.size()), 2*n);
        for (uint j=1; j<=niters; j++) {
            power = power * b;
            result.add_(power * Complex(pow(-1, j+1) / j, 0.));
        }
        return result;
    }
};
#endif