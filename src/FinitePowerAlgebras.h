#ifndef ALGEBRAS_H
#define ALGEBRAS_H
#include "BaseAlgebraRelations.h"
using namespace std; 

// Exterior and Clifford Algebras

// General design for a particular algebra:
    // Define a commutation relation accepting conjugation relation as template argument 
    // Instantiate a template alias

// Exterior algebra: 
// Basic commutation relation for exterior algebra: 
//    All generators anti-commute, and >1 powers annihilate
//    Accepts variable conjugation relation
template<typename ConjRelation>
struct ExteriorCommRelation: public ConjRelation {
    CoeffMap commute_noncanonical(uint i, uint j) const override {
        CoeffMap result;
        result[{j, i}] = Complex(-1., 0.);
        return result; 
    }
    
    tuple<uint, Complex> homogeneous_exponent(uint gidx, uint pow) const override {
        static_cast<void>(gidx);
        // Annihilates powers > 1
        auto scale = pow > 1 ? Complex(0., 0.) : Complex(1., 0.);
        return tuple<uint, Complex>(pow, scale);
    }
};


// Basic operations for exterior algebra:
    // Logarithms and exponentials defined via power series
template<typename AlgRelation>
class ExteriorAlgebra : public BaseAlgebra<AlgRelation> {
public:
    using Element = AlgebraElement<AlgRelation>;
    // Grassmann exponentiation: approximation
    Element exp(const Element& a) {
        auto power = a.one(), result=a.one(); // Start with the multiplicative identity;
        uint niters = 20;
        for (uint j=1; j<=niters; j++) {
            power = power * a;
            result.add_(power * Complex(1./static_cast<double>(factorial(j)), 0));
        }
        return result;
    }

    // Taylor expansion definition of the logarithm. 
    Element log(const Element& a) {
        auto b = a - Complex(1., 0.);
        if (b.coeffs.contains(KeyType(AlgRelation::num_generators(), 0u))) {
            throw(std::invalid_argument(
                "Exterior logarithm of " + b.to_string() 
                + " not guaranteed to converge"));
        }
        auto power = a.one(), result = a.zero();
        uint niters = min(static_cast<uint>(b.coeffs.size()), AlgRelation::num_generators());
        for (uint j=1; j<=niters; j++) {
            power = power * b;
            result.add_(power * Complex(pow(-1, j+1) / j, 0.));
        }
        return result;
    }
};


template<uint n>
using ExtFreeConjAlgebra = ExteriorAlgebra<ExteriorCommRelation<FreeConjRelation<n>>>;
template<uint n>
using ExtSelfConjAlgebra = ExteriorAlgebra<ExteriorCommRelation<SelfConjRelation<n>>>;


// Clifford algebra: anticommute with each other and square to one
//  conjugates to self
template<typename ConjRelation>
struct CliffordCommRelation: public ConjRelation {
    CoeffMap commute_noncanonical(uint i, uint j) const override {
        CoeffMap result;
        result[{j, i}] = Complex(-1., 0.);
        return result; 
    }
    tuple<uint, Complex> homogeneous_exponent(uint gidx, uint pow) const override {
        static_cast<void>(gidx);
        return tuple<uint, Complex>(pow % 2, Complex(1., 0.));
    }
};

template<uint n>
using CliffordAlgebra = BaseAlgebra<CliffordCommRelation<SelfConjRelation<n>>>;


// CAR algebra using fermionic creation / annihilation operators
// Canonical ordering {a1, a1*, a2, a2*, ...}
//    The only nontrivially anti-commuting relation is a1* a1 = 1 - a1 a1* (i=1, j=0)
template<uint n>
struct CanonicalAnticommRelation: public FreeConjRelation<n> {
    CoeffMap commute_noncanonical(uint i, uint j) const override {
        CoeffMap result;
        if (i%2==1 && j==i-1) {
            result[{}] = Complex(1., 0.);
        }
        result[{j, i}] = Complex(-1., 0.);
        return result; 
    }
    tuple<uint, Complex> homogeneous_exponent(uint gidx, uint pow) const override {
        static_cast<void>(gidx);
        // Annihilates powers > 1
        auto scale = pow > 1 ? Complex(0., 0.) : Complex(1., 0.);
        return tuple<uint, Complex>(pow, scale);
    }
};

template<uint n> 
using CARAlgebra = BaseAlgebra<CanonicalAnticommRelation<n>>;

#endif