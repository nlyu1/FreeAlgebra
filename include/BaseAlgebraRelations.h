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
    // Specifies the trace of generators. 
    //    Traces are expected to extend multiplicatively across different generators
    virtual Complex tr(uint gidx, uint pow) const {
        static_cast<void>(gidx);
        static_cast<void>(pow);
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
    // Given the algebra-specific generator trace relation, 
    //    computes the trace of a monomial
    virtual Complex monomial_tr(const KeyType& monomial) {
        Complex result = Complex(0.);
        for (uint j=0; j<monomial.size(); j++) {
            if (result == Complex(0., 0.)) {
                return result; 
            }
            result *= tr(j, monomial[j]); 
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


// Assumes that 1 is represented as id on d-dimensional 
//    space and all other operators are traceless. 
template<uint n, uint d>
struct ScalarTraceRelation: virtual public BaseRelation<n> {
    // Specifies the trace of generators  
    Complex tr(uint gidx, uint pow) const override {
        static_cast<void>(gidx);
        return pow == 0 ? Complex(static_cast<double>(d)) : Complex(0., 0.);
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
            CoeffMap({{KeyType(Relation::num_generators(), 0u), Complex(1., 0.)}})
        );
    }

    Element commutator(const Element& a, const Element& b){
        return a * b - b * a;
    }

    Element anticommutator(const Element& a, const Element& b){
        return a * b + b * a;
    }
};




// Complex template parameter
// Right now double-types are not supported, so we use int divided by 1e7
template<int real, int imag>
struct ComplexParameter {
    static constexpr Complex value{real/10000000, imag/10000000};
};
// Pass these into the third template argument of ProdAlgebra 
//    to obtain commuting (resp. anticommuting) tensor product algebras
constexpr ComplexParameter<10000000, 0> PROD_COMMUTE;
constexpr ComplexParameter<-10000000, 0> PROD_ANTICOMMUTE;


// Takes the generators of two constituent relations commute with each other 
//    up to a complex phase specified by CommutePhase
template<typename LRel, typename RRel, typename CommutePhase>
struct ProductRelation: 
    public BaseRelation<LRel::num_generators() + RRel::num_generators()> { 
    static LRel& Lrel() {
        static LRel Lrel_{};
        return Lrel_;
    }
    static RRel& Rrel() {
        static RRel Rrel_{};
        return Rrel_;
    }

    // The return-type should be power representation 
    CoeffMap commute_noncanonical(uint i, uint j) const override {
        auto pivot = LRel::num_generators();
        CoeffMap result; 
        if (j < pivot && i >= pivot) { // j < pivot <= i
            /// Replace this Complex(-1., 0.) with the third template argument
            result[{j, i}] = CommutePhase::value; 
        } else if (i < pivot) { 
            // then j < i < pivot
            //  Compute CR using RRel then extend to the whole algebra
            result = Lrel().commute_noncanonical(i, j);
        } else if (j >= pivot) { 
            // pivot <= j < i 
            //   Compute CR using LRel then extend to the whole algebra
            auto result_ = Rrel().commute_noncanonical(i-pivot, j-pivot); 
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
        auto pivot = LRel::num_generators();
        return (i<pivot) ? Lrel().conj(i) : Rrel().conj(i - pivot) + pivot;
    }

    tuple<uint, Complex> homogeneous_exponent(uint gidx, uint pow) const override {
        auto pivot = LRel::num_generators();
        if (gidx < pivot) {
            return Lrel().homogeneous_exponent(gidx, pow);
        } else {
            return Rrel().homogeneous_exponent(gidx - pivot, pow);
        }
    }

    Complex tr(uint gidx, uint pow) const override {
        return gidx < LRel::num_generators()? 
            Lrel().tr(gidx, pow) : Rrel().tr(gidx - LRel::num_generators(), pow);
    }
};


template<typename LAlg, typename RAlg, typename CommutePhase>
class ProductAlgebra: 
    public BaseAlgebra<ProductRelation<
        typename LAlg::Relation, typename RAlg::Relation, CommutePhase
    >> {
public:
    using LRel = typename LAlg::Relation; 
    using RRel = typename RAlg::Relation;

    using LElm = typename LAlg::Element; 
    using RElm = typename RAlg::Element; 
    using Relation = ProductRelation<typename LAlg::Relation, 
                    typename RAlg::Relation, CommutePhase>; 
    using Element = typename BaseAlgebra<Relation>::Element; 

    // Takes the tensor product of x and y
    Element kron(const LElm& x, const LElm& y) {
        CoeffMap coeffs; 
        for (auto const& pair1:x.coeffs) {
            for (auto const& pair2:y.coeffs) {
                coeffs[mergeVectors(
                    pair1.first, pair2.first)] = pair1.second * pair2.second; 
            }
        }
        return Element(coeffs);
    }

    // Extend an element of the left algebra to the right 
    Element extR(const LElm& x) {
        CoeffMap coeffs;
        for (auto const& pair:x.coeffs) {
            coeffs[mergeVectors(
                    pair.first, 
                    KeyType(RRel::num_generators(), 0u)
                )] = pair.second; 
        }
        return Element(coeffs);
    }

    // Extend an element of the right algebra to the left
    Element extL(const RElm& x) {
        CoeffMap coeffs;
        for (auto const& pair:x.coeffs) {
            coeffs[mergeVectors(
                    KeyType(LRel::num_generators(), 0u),
                    pair.first
                )] = pair.second; 
        }
        return Element(coeffs);
    }

    // Project a product element to the left algebra. 
    //    Any component with nontrivial support on the right altegra is thrown away
    LElm projL(const Element& x) {
        CoeffMap coeffs;
        for (auto const& pair:x.coeffs) {
            // The new key is the left "half" of the monomial
            auto pivot = pair.first.begin()+LRel::num_generators();
            auto keyL = KeyType(pair.first.begin(), pivot);
            auto keyR = KeyType(pivot, pair.first.end());
            coeffs[keyL] = pair.second * static_cast<double>(allzero(keyR));
        }
        return LElm(coeffs);
    }

    RElm projR(const Element& x) {
        CoeffMap coeffs;
        for (auto const& pair:x.coeffs) {
            auto pivot = pair.first.begin()+LRel::num_generators();
            auto keyL = KeyType(pair.first.begin(), pivot);
            auto keyR = KeyType(pivot, pair.first.end());
            coeffs[keyR] = pair.second * allzero(keyL);
        }
        return RElm(coeffs);
    }

    // Partial trace over the right algebra
    LElm trR(const Element& x) {
        CoeffMap coeffs;
        for (auto const& pair:x.coeffs) {
            // The new key is the left "half" of the monomial
            auto pivot = pair.first.begin()+LRel::num_generators();
            auto keyL = KeyType(pair.first.begin(), pivot);
            auto keyR = KeyType(pivot, pair.first.end());
            coeffs[keyL] = pair.second * Relation::Rrel().monomial_tr(keyR);
        }
        return LElm(coeffs);
    }

    RElm trL(const Element& x) {
        CoeffMap coeffs;
        for (auto const& pair:x.coeffs) {
            auto pivot = pair.first.begin()+LRel::num_generators();
            auto keyL = KeyType(pair.first.begin(), pivot);
            auto keyR = KeyType(pivot, pair.first.end());
            coeffs[keyR] = pair.second * Relation::Lrel().monomial_tr(keyL);
        }
        return RElm(coeffs);
    }
};
#endif