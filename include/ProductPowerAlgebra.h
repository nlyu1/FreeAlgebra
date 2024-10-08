#ifndef POWER_H
#define POWER_H
#include "Elements.h"
#include "BaseAlgebraRelations.h"
using namespace std;

// FieldType template parameter
// Right now double-types are not supported, so we use int divided by 1e7
template<int real, int imag>
struct ComplexParameter {
    static constexpr std::complex<double> value{real/10000000., imag/10000000.};
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
            /// Replace this FieldType(-1., 0.) with the third template argument
            result[{j, i}] = FieldType(CommutePhase::value.real(), CommutePhase::value.imag()); 
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

    FieldType tr(uint gidx, uint pow) const override {
        return gidx < LRel::num_generators()? 
            Lrel().tr(gidx, pow) : Rrel().tr(gidx - LRel::num_generators(), pow);
    }

    std::string to_string(uint v) const override {
        if (v < LRel::num_generators()) {
            return Lrel().to_string(v) + std::string("L"); 
        } else {
            return Lrel().to_string(v - LRel::num_generators())
                + std::string("R"); 
        }
    }
};


template<typename LAlg, typename RAlg, typename CommutePhase>
class ProductAlgebra: 
    public BaseAlgebra<ProductRelation<
        typename LAlg::Relation, typename RAlg::Relation, CommutePhase>
    > {
public:
    using LRel = typename LAlg::Relation; 
    using RRel = typename RAlg::Relation;

    using LElm = typename LAlg::Element; 
    using RElm = typename RAlg::Element; 
    using LAlgebra = LAlg; 
    using RAlgebra = RAlg; 
    using Relation = ProductRelation<typename LAlg::Relation, 
                    typename RAlg::Relation, CommutePhase>; 
    using Element = typename BaseAlgebra<Relation>::Element; 

    static LAlgebra& lAlg() {
        static LAlgebra lAlg_{};
        return lAlg_;
    }

    static RAlgebra& rAlg() {
        static RAlgebra rAlg_{};
        return rAlg_;
    }

    // Takes the tensor product of x and y
    Element kron(const LElm& x, const RElm& y) {
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
            auto scalar = pair.second * Relation::Rrel().monomial_tr(keyR);
            if (scalar == FieldType(0.)) {
                continue; 
            }
            // cout << prettyPrint(keyL) << " " << prettyPrint(keyR) << ": " 
            // << pair.second << ", " << Relation::Rrel().monomial_tr(keyR)
            // << endl;
            coeffs[keyL] = scalar;
        }
        auto result = LElm(coeffs);
        cout << result << endl;
        return result;
    }

    RElm trL(const Element& x) {
        CoeffMap coeffs;
        for (auto const& pair:x.coeffs) {
            auto pivot = pair.first.begin()+LRel::num_generators();
            auto keyL = KeyType(pair.first.begin(), pivot);
            auto keyR = KeyType(pivot, pair.first.end());
            auto scalar = pair.second * Relation::Lrel().monomial_tr(keyL);
            if (scalar == FieldType(0.)) {
                continue; 
            }
            // cout << prettyPrint(keyL) << " " << prettyPrint(keyR) << ": " 
            // << pair.second << ", " << Relation::Lrel().monomial_tr(keyL)
            // << endl;
            coeffs[keyR] = pair.second * Relation::Lrel().monomial_tr(keyL);
        }
        return RElm(coeffs);
    }
};


typedef std::vector<FieldType> ComplexDispType;
template<typename T>
std::vector<T> add_vec(const std::vector<T>& vec1, 
    const std::vector<T>& vec2) {
    size_t len = std::min(vec1.size(), vec2.size());
    std::vector<T> result;
    for (size_t i = 0; i < len; ++i) {
        result.push_back(vec1[i] + vec2[i]);
    }
    return result;
}
template<typename T>
std::vector<T> neg_vec(const std::vector<T>& v) {
    std::vector<T> result;
    for (auto x:v) {
        result.push_back(x * -1);
    }
    return result;
}


template<typename LAlg, typename RAlg>
class ExtProductAlgebra: 
    public ProductAlgebra<RAlg, RAlg, decltype(PROD_ANTICOMMUTE)> {
public:
    using LRel = typename LAlg::Relation; 
    using RRel = typename RAlg::Relation;

    using LElm = typename LAlg::Element; 
    using RElm = typename RAlg::Element; 
    using LAlgebra = LAlg; 
    using RAlgebra = RAlg; 
    using Relation = ProductRelation<typename LAlg::Relation, 
                    typename RAlg::Relation, decltype(PROD_ANTICOMMUTE)>; 
    using Element = typename BaseAlgebra<Relation>::Element; 

    static LAlgebra& lAlg() {
        static LAlgebra lAlg_{};
        return lAlg_;
    }

    static RAlgebra& rAlg() {
        static RAlgebra rAlg_{};
        return rAlg_;
    }

    LElm dR(const Element& x, 
        const ComplexDispType& c=ComplexDispType(LAlgebra::num_generators(), FieldType(1.))) {
        CoeffMap coeffs; 
        assert (c.size() == RAlgebra::num_generators()); 
        auto s = FieldType(1.);
        for (auto v:c) {
            s = s * v.absq(); 
        }
        for (auto const& p:x.coeffs) {
            // Extract the left half of the coefficients 
            auto scale = p.second; 
            auto LKeys = KeyType(p.first.begin(), 
                p.first.begin()+LAlgebra::num_generators());
            auto RKeys = KeyType(p.first.begin()+LAlgebra::num_generators(), 
                p.first.end());
            scale = scale * rAlg().td(RElm({{RKeys, FieldType(1.)}})); 
            if (scale == 0.) continue; 
            coeffs[LKeys] = scale / s; 
        }
        return LElm(coeffs);
    }
};

#endif