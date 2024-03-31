#ifndef ALGEBRAS_H
#define ALGEBRAS_H
#include "element.h"
using namespace std; 

// Algebraic relations for exterior algebra elements with conjugation
// Canonical ordering (a1, a1*, a2, a2*, ...)
// Everything about the algebra is encoded in the following struct 
//    templated on the number of generators 
template<uint n>
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

    static constexpr uint num_generators() {
        return n;
    }

    static std::string to_string(uint v) {
        return v%2==0 ? std::to_string(v/2) : std::to_string(v/2)+"*";
    }
};

// Dirac commutation relation (CAR)
// Canonical ordering (a1, a1*, a2, a2*, ...)
template<uint n>
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

    static constexpr uint num_generators() {
        return n;
    }
};


// The algebra itself is templated on 
//   the number of modes (half the number of generators)
template<uint n>
class ExtConjAlg {
public:
    ExtConjAlg() = default;

    // Returns the corresponding annihilator
    AlgebraElement<ExtConjRelation<2*n>> operator()(uint i) const {
        assert(i >= 0 && i<n); 
        KeyType coeffs({2*i});
        return AlgebraElement<ExtConjRelation<2*n>>({{gen_to_power_repr(coeffs, 2*n), 1}});
    }

    // Grassmann exponentiation: only need to compute a few terms, yeah
    AlgebraElement<ExtConjRelation<2*n>> exp(const AlgebraElement<ExtConjRelation<2*n>>& a) {
        auto power = a.one(), result=a.one(); // Start with the multiplicative identity;
        uint niters = min(static_cast<uint>(a.coeffs.size()), 2*n);
        for (uint j=1; j<=niters; j++) {
            power = power * a;
            result.add_(power * Complex(1./static_cast<double>(factorial(j)), 0));
        }
        return result;
    }

    // Taylor expansion definition of the logarithm. 
    AlgebraElement<ExtConjRelation<2*n>> log(const AlgebraElement<ExtConjRelation<2*n>>& a) {
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