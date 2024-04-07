#ifndef ALGEBRA_TEMPL_H
#define ALGEBRA_TEMPL_H
#include "Elements.h"
using namespace std;
// Generic algebra constructs: 
//      Relations (BaseRelation and common conjugate relations)
//      Algebras (BaseAlgebra and Alternating product relations)

// Most basic relation: defines the num_generators() and to_string
template<uint n>
struct BaseRelation {
    /// Algebra-specific relations
    // Specifies commutation relation when i>j. Algebra-specific
    virtual CoeffMap commute_noncanonical(uint i, uint j) const {
        static_cast<void>(i); static_cast<void>(j);
        throw std::invalid_argument("Noncanonical commutation relation needs to be defined");
    }
    // Specifies conjugation relation on a generator 
    virtual uint conj(uint i) const {
        static_cast<void>(i);
        throw std::invalid_argument("Conjugation not implemented");
    }
    // Specifies the trace of generators. 
    //    Traces are expected to extend multiplicatively across different generators
    virtual FieldType tr(uint gidx, uint pow) const {
        static_cast<void>(gidx);
        static_cast<void>(pow);
        throw std::invalid_argument("Tracial relation not implemented");
    }

    /// Default implementations 
    CoeffMap commute(uint i, uint j) const {
        assert (i>=0 && i<n);
        assert (j>=0 && j<n);
        if (i < j) {
            throw std::invalid_argument("Expects non-canonical ordering");
        }
        return commute_noncanonical(i, j);
    }
    virtual std::string to_string(uint v) const {
        return std::to_string(v);
    }
    static constexpr uint num_generators() {
        return n;
    }
    // Given the algebra-specific generator trace relation, 
    //    computes the trace of a monomial
    virtual FieldType monomial_tr(const KeyType& monomial) {
        FieldType result = FieldType(1.);
        for (uint j=0; j<monomial.size(); j++) {
            if (result == FieldType(0., 0.)) {
                return result; 
            }
            result = result * tr(j, monomial[j]); 
        }
        return result; 
    }
    virtual ~BaseRelation() = default; // A warning technicality
};


// Free-conjugation relation: formal conjugates are 
// Canonical ordering (a1, a1*, a2, a2*, ..., a(n/2), a(n/2)*)
// n must be even 
template<uint n>
struct FreeConjRelation: virtual public BaseRelation<n> {
    static_assert(n % 2 == 0, 
        "Algebra with formal conjugate generators contain need an even number of generators");
    // Specifies which element to conjugate to
    uint conj(uint i) const override {
        return (i % 2 == 0) ? i+1 : i-1;
    }
    std::string to_string(uint v) const override {
        return std::to_string(v / 2) + (v % 2 == 0 ? "" : "*");
    }
};

#define TRACE_SCALE 1000000000

// Assumes that 1 is represented as id on d-dimensional 
//    space and all other operators are traceless. 
template<uint n, uint d>
struct ScalarTraceRelation: virtual public BaseRelation<n> {
    // Specifies the trace of generators  
    FieldType tr(uint gidx, uint pow) const override {
        static_cast<void>(gidx);
        return pow == 0 ? FieldType(static_cast<double>(d) / TRACE_SCALE) : FieldType(0., 0.);
    }
};

// Basic self-conjugate relation: each generator conjugates to itself
// Canonical ordering (a1, a2, ...)
template<uint n>
struct SelfConjRelation: virtual public BaseRelation<n> {
    // Specifies which element to conjugate to
    uint conj(uint i) const override {
        return i;
    }
    std::string to_string(uint v) const override {
        return std::to_string(v);
    }
};

template<typename AlgRelation>
class BaseAlgebra {
public:
    using Element = AlgebraElement<AlgRelation>; 
    using Relation = AlgRelation; 
    
    BaseAlgebra<AlgRelation>() = default;

    // Returns the corresponding annihilator
    Element operator()(uint i) const {
        assert(i >= 0 && i< Relation::num_generators()); 
        KeyType coeffs({i});
        return Element({{gen_to_power_repr(coeffs, Relation::num_generators()), 1}});
    }

    Element zero() const {
        return Element(CoeffMap({}));
    }

    Element one() const {
        return Element(
            CoeffMap({{KeyType(Relation::num_generators(), 0u), FieldType(1., 0.)}})
        );
    }

    Element commutator(const Element& a, const Element& b){
        return a * b - b * a;
    }

    Element anticommutator(const Element& a, const Element& b){
        return a * b + b * a;
    }

    static constexpr uint num_generators() {
        return Relation::num_generators();
    }
};

#endif