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
template<uint n> 
struct ExteriorCommRelation: virtual public BaseRelation<n> {
    CoeffMap commute_noncanonical(uint i, uint j) const override {
        CoeffMap result;
        // Squares annihilate
        if (i == j) {
            return result; 
        }
        result[{j, i}] = FieldType(-1., 0.);
        return result; 
    }
};



template<uint n>
struct ExtFreeConjRelation: 
    public ExteriorCommRelation<n>, public FreeConjRelation<n>, public ScalarTraceRelation<n, 2>
    {};
template<uint n>
struct ExtSelfConjRelation: 
    public ExteriorCommRelation<n>, public SelfConjRelation<n>, public ScalarTraceRelation<n, 2>
    {};

// Uses exterior commutation relation and free conjugation relation 
template<uint n>
using ExtFreeConjAlgebra = BaseAlgebra<ExtFreeConjRelation<n>>;
template<uint n>
using ExtSelfConjAlgebra = BaseAlgebra<ExtSelfConjRelation<n>>;


// Clifford algebra: anticommute with each other and square to one
//  conjugates to self
template<uint n>
struct CliffordCommRelation: virtual public BaseRelation<n> {
    CoeffMap commute_noncanonical(uint i, uint j) const override {
        CoeffMap result;
        if (i == j) {
            result[{}] = FieldType(1.);
            return result;
        }
        result[{j, i}] = FieldType(-1.);
        return result; 
    }
};

template<uint n>
struct CliffordRelation:
    public CliffordCommRelation<n>, public SelfConjRelation<n>, public ScalarTraceRelation<n, 2>
    {};

template<uint n>
using CliffordAlgebra = BaseAlgebra<CliffordRelation<n>>;


// CAR algebra using fermionic creation / annihilation operators
// Canonical ordering {a1, a1*, a2, a2*, ...}
//    The only nontrivially anti-commuting relation is a1* a1 = 1 - a1 a1* (i=1, j=0)
template<uint n>
struct CanonicalAnticommRelation: public FreeConjRelation<n> {
    CoeffMap commute_noncanonical(uint i, uint j) const override {
        CoeffMap result;
        if (i == j) {
            return result; 
        }
        if (i%2==1 && j==i-1) {
            result[{}] = FieldType(1., 0.);
        }
        result[{j, i}] = FieldType(-1., 0.);
        return result; 
    }
};

template<uint n> 
using CARAlgebra = BaseAlgebra<CanonicalAnticommRelation<n>>;

#endif