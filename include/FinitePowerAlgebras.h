#ifndef ALGEBRAS_H
#define ALGEBRAS_H
#include "BaseAlgebraRelations.h"
using namespace std; 

static const uint sqrt_two = static_cast<uint>(TRACE_SCALE * 1.4142135623730951);
static const uint two = static_cast<uint>(TRACE_SCALE * 2);

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


// The trace of exterior commutation is not well-defined 
//    by itself since its generators are not tensor-independent
template<uint n>
struct ExtFreeConjRelation: 
    public ExteriorCommRelation<n>, public FreeConjRelation<n>
    {};
template<uint n>
struct ExtSelfConjRelation: 
    public ExteriorCommRelation<n>, public SelfConjRelation<n>
    {};

// Augments exterior differentiation 
template<typename AlgRelation>
class ExteriorAlgebra: public BaseAlgebra<AlgRelation> {
public:
    using BaseAlgebra = BaseAlgebra<AlgRelation>;
    using Element = AlgebraElement<AlgRelation>; 
    using Relation = AlgRelation; 

    static Relation& rel() {
        static Relation rel_{};
        return rel_;
    }

    static BaseAlgebra& alg() {
        static BaseAlgebra alg_{};
        return alg_;
    }

    // Derivative operator (equivalent to integral operator)
    Element d(const Element& x, uint idx) const {
        CoeffMap coeffs; 
        for (auto const& p: x.coeffs) {
            // Only consider the term if it contains the required variable 
            if (p.first[idx] == 0) continue; 
            FieldType scale = p.second;
            // For each preceding generator, add a negative sign 
            for (uint i=0; i<idx; i++) {
                if (p.first[i] != 0) scale = scale * -1.;
            } 
            KeyType k(p.first);
            k[idx] = 0;
            coeffs[k] = scale; 
        }
        return Element(coeffs);
    }

    // Total derivative: apply dx(0) then dx(1) then ...
    FieldType td(const Element& x) const {
        auto result = x;
        for (uint i=0; i<Relation::num_generators(); i++) {
            result = d(result, i);
        }
        if (result.coeffs.size() == 0) {
            return FieldType(0.);
        } else if (result.coeffs.size() == 1) {
            return result.coeffs.begin() -> second;
        } else {
            throw(std::invalid_argument("Unexpected case during total derivative"));
        }
    }

    // Logarithms are always guaranteed to exist
    // Shift the center so that it converges over a finite number of terms 
    Element log(const Element& x) {
        // The key for scalar offshift 
        auto onekey = KeyType(Relation::num_generators(), 0u);
        double center = 0.;
        if (x.coeffs.contains(onekey)) {
            auto c = x.coeffs.at(onekey);
            if (std::abs(c.imag()) > 1e-10 || std::abs(c.real()) < 1e-10) {
                throw(std::invalid_argument(fmt::format(
                    "Expected vanishing imaginary offset and nontrivial real offset but received {}\n", 
                    c.to_string()
                )));
            }
            center = c.real(); 
        } else {
            throw(std::invalid_argument(
                    "Expected nontrivial real offset for exterior logarithm.\n"
                ));
        }
        // Subtract the center so only grassmann elements exist 
        auto b = x / center - 1; 
        auto power = alg().one(), result=alg().zero()+std::log(center);
        for (uint j=1; j<=Relation::num_generators(); j++) {
            // cout << j << endl;
            power = power * b;
            if (power == 0.) break;
            result.add_(power*(std::pow(-1, j+1)/j));
        }
        // Consistency check for logarithm
        cout << "Log-difference: " << (x - result.exp()).norm() << endl;
        // assert ((x - result.exp()).norm() < 1e-10); 
        return result; 
    }
};


template<uint n>
using ExtFreeConjAlgebra = ExteriorAlgebra<ExtFreeConjRelation<n>>;

template<uint n>
using ExtSelfConjAlgebra = ExteriorAlgebra<ExtSelfConjRelation<n>>;

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
    public CliffordCommRelation<n>, public SelfConjRelation<n>, 
    public ScalarTraceRelation<n, sqrt_two>
    {};

template<uint n>
using CliffordAlgebra = BaseAlgebra<CliffordRelation<n>>;


// CAR algebra using fermionic creation / annihilation operators
// Canonical ordering {a1, a1*, a2, a2*, ...}
//    The only nontrivially anti-commuting relation is a1* a1 = 1 - a1 a1* (i=1, j=0)
template<uint n>
struct CanonicalAnticommRelation: 
    public FreeConjRelation<n>, public ScalarTraceRelation<n, TRACE_SCALE> {
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