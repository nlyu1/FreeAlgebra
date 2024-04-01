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
    // Specifies (gidx)**pow = (Complex) * (gidx)**(uint)
    virtual tuple<uint, Complex> homogeneous_exponent(uint gidx, uint pow) const {
        static_cast<void>(gidx); static_cast<void>(pow);
        throw std::invalid_argument("Noncanonical commutation relation needs to be defined");
    }
    // Specifies conjugation relation on a generator 
    virtual uint conj(uint i) const {
        static_cast<void>(i);
        throw std::invalid_argument("Conjugation not implemented");
    }

    /// Default implementations 
    CoeffMap commute(uint i, uint j) const {
        assert (i>=0 && i<n);
        assert (j>=0 && j<n);
        if (i <= j) {
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
    virtual ~BaseRelation() = default; // A warning technicality
};


// Free-conjugation relation: formal conjugates are 
// Canonical ordering (a1, a1*, a2, a2*, ..., a(n/2), a(n/2)*)
// n must be even 
template<uint n>
struct FreeConjRelation: public BaseRelation<n> {
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

// Basic self-conjugate relation: each generator conjugates to itself
// Canonical ordering (a1, a2, ...)
template<uint n>
struct SelfConjRelation: public BaseRelation<n> {
    // Specifies which element to conjugate to
    uint conj(uint i) const override {
        return i;
    }
    std::string to_string(uint v) const override {
        return std::to_string(v);
    }
};

// Takes the generators of two constituent relations to be opposite each other
template<typename RelA, typename RelB>
struct AltProductRelation : 
    public BaseRelation<RelA::num_generators() + RelB::num_generators()> { 

    // The return-type should be power representation 
    CoeffMap commute_noncanonical(uint i, uint j) const override {
        auto pivot = RelA::num_generators();
        CoeffMap result; 
        if (j < pivot && i >= pivot) { // j < pivot <= i
            result[{j, i}] = Complex(-1., 0.);
        } else if (i < pivot) { 
            // then j < i < pivot
            //  Compute CR using RelB then extend to the whole algebra
            result = RelA().commute_noncanonical(i, j);
        } else if (j >= pivot) { 
            // pivot <= j < i 
            //   Compute CR using RelA then extend to the whole algebra
            auto result_ = RelB().commute_noncanonical(i-pivot, j-pivot); 
            // For each entry in the original cr result
            for (auto& pair : result_) {
                KeyType key_(pair.first);
                for (size_t k=0; k<pair.first.size(); k++) {
                    key_[k] += pivot;
                }
                result[key_] = pair.second;
            }
        } else {
            throw(std::invalid_argument("Product relation unexpected"));
        }
        return result; 
    }

    uint conj(uint i) const override {
        auto pivot = RelA::num_generators();
        return (i<pivot) ? RelA().conj(i) : RelB().conj(i - pivot) + pivot;
    }

    tuple<uint, Complex> homogeneous_exponent(uint gidx, uint pow) const override {
        auto pivot = RelA::num_generators();
        if (gidx < pivot) {
            return RelA().homogeneous_exponent(gidx, pow);
        } else {
            return RelB().homogeneous_exponent(gidx - pivot, pow);
        }
    }
};


template<typename AlgRelation>
class BaseAlgebra {
public:
    using Element = AlgebraElement<AlgRelation>; 
    BaseAlgebra<AlgRelation>() = default;

    // Returns the corresponding annihilator
    Element operator()(uint i) const {
        assert(i >= 0 && i< AlgRelation::num_generators()); 
        KeyType coeffs({i});
        return Element({{gen_to_power_repr(coeffs, AlgRelation::num_generators()), 1}});
    }

    Element zero() const {
        return Element(CoeffMap({}));
    }

    Element one() const {
        return Element(
            CoeffMap({{KeyType(AlgRelation::num_generators(), 0u), Complex(1., 0.)}})
        );
    }

    Element commutator(const Element& a, const Element& b){
        return a * b - b * a;
    }

    Element anticommutator(const Element& a, const Element& b){
        return a * b + b * a;
    }
};
#endif