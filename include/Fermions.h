#ifndef FERMIONS_H
#define FERMIONS_H
#include <cmath>
#include <algorithm>
#include "BaseAlgebraRelations.h"
#include "ProductPowerAlgebra.h"
// #include "FinitePowerAlgebras.h"
#include "utils.h"
using namespace std; 

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
// ComplexDispType symplectic_product(const ComplexDispType& a, const ComplexDispType& b) {
//     FieldType result(0.);
//     assert (a.size() == b.size());
    
//     for (uint i=0; i<a.size(); i++) {
//         result = result + 
//     }
// }

template<uint n>
class DiracAlgebra: public ProductAlgebra<
        CARAlgebra<2*n>, ExtFreeConjAlgebra<2*n>, decltype(PROD_ANTICOMMUTE)> {
public:
    using BaseAlg = ProductAlgebra<
        CARAlgebra<2*n>, ExtFreeConjAlgebra<2*n>, decltype(PROD_ANTICOMMUTE)>; 
    using Element = BaseAlg::Element; 

    using LAlgebra = BaseAlg::LAlgebra; 
    using RAlgebra =  BaseAlg::RAlgebra; 
    using LElm = typename LAlgebra::Element; 
    using RElm = typename RAlgebra::Element; 
    // The majorana algebra, in another basis
    using MAlgebra = CliffordAlgebra<2*n>;
    using MElm = typename MAlgebra::Element;

    static BaseAlg& alg() {
        static BaseAlg alg_{};
        return alg_;
    }

    static LAlgebra& lAlg() {
        static LAlgebra lAlg_{};
        return lAlg_;
    }

    static RAlgebra& rAlg() {
        static RAlgebra rAlg_{};
        return rAlg_;
    }

    static MAlgebra& mAlg() {
        static MAlgebra mAlg_{};
        return mAlg_;
    }

    // Defines the vacuum state 
    Element vac() const {
        auto result = alg().one();
        for (uint i=0; i<n; i++) {
            result = result * (a(i) * a(i).conj());
        }
        return result; 
    }

    // The extended annihilation operators 
    Element a(uint i) const {
        assert (i>=0 && i<n); 
        return alg().extR(lAlg()(2*i));
    }

    // Returns the displacement operator 
    Element disp(ComplexDispType d) const {
        assert (d.size() == n);
        auto result = alg().zero();
        auto dvec = disp_vec(d); 
        for (uint i=0; i<n; i++) {
            result.add_(a(i).conj() * dvec[i] + a(i) * dvec[i].conj()); 
        }
        return result.exp(); 
    }

    // Encoding of a double element as anticommuting exterior vector
    std::vector<Element> disp_vec(const ComplexDispType& a) const {
        auto r = rAlg();
        assert (a.size() == n);
        std::vector<Element> ans; 
        for (uint i=0; i<n; i++) {
            ans.push_back(alg().extL(r(2*i) * a[i]));
        }
        return ans; 
    }

    // Commutator corresponding to displacement addition
    // a.disp(add_vec(d1, d2)) == D1 * D2 * a.disp_add(d1, d2)
    Element disp_add(const ComplexDispType& a, const ComplexDispType& b) const {
        auto avec = disp_vec(a), bvec = disp_vec(b); 
        auto k = alg().zero();
        for (uint i=0; i<n; i++) {
            k.add_(avec[i] * bvec[i].conj() + avec[i].conj() * bvec[i]);
        }
        return (k / 2.).exp();
    }

    Element coherent(const ComplexDispType& a) {
        auto D = disp(a); 
        return D * vac() * D.conj(); 
    }

    // Integration against a vector of Grassmann values
    //   The grassmann values themselves are by default with unity coefficients
    LElm dR(const Element& x, 
        const ComplexDispType& c=ComplexDispType(n, FieldType(1.))) {
        CoeffMap coeffs; 
        assert (c.size() == n); 
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

    // Given real coefficients for majorana generators, 
    //    ensures that the result is Hermitian and assembles the Hamiltonian 
    Element hamiltonian(const RealCoeffMap& coeffs_) const {
        CoeffMap coeffs;
        for (const auto& pair: coeffs_) {
            auto degree = pair.first.size();
            auto key = gen_to_power_repr(pair.first, 2*n);
            // Check that we don't have the wrong representation
            for (auto const k: key) {
                if (k>1) {
                    throw(std::invalid_argument(fmt::format(
                        "Expected mult-representation but obtained power representation {}", 
                        prettyPrint(pair.first)
                    )));
                }
            }
            if (degree * (degree - 1) / 2 % 2 == 0){ 
                coeffs[key] = pair.second * M_PI;
            } else {
                coeffs[key] = FieldType(0, 1.) * pair.second * M_PI;
            }
        }
        auto majorana_result = MElm(coeffs);
        // assert (majorana_result == car_to_majorana(majorana_to_car(majorana_result)));
        auto result = alg().extR(majorana_to_car(majorana_result));
        assert (result.isreal()); 
        return result;
    }

    // Assembles a unitary from the coefficient of its exponents 
    Element unitary(const RealCoeffMap& coeffs) const {
        return (hamiltonian(coeffs) * FieldType(0., 1.)).exp();
    }

    // Assembles a thermal state from its Hamiltonian 
    Element density(const RealCoeffMap& coeffs) const {
        auto expH = hamiltonian(coeffs).exp();
        return expH / tr(expH).coeffs.begin()->second;
    }

    // Uses coeffs to compute a unitary acting on the zero state
    Element pure(const RealCoeffMap& coeffs) const {
        auto U = unitary(coeffs); 
        return U * vac() * U.conj();
    }

    // Uses coeffs to compute a unitary acting on the zero state
    Element pure_gaussian(const RealCoeffMap& coeffs) const {
        for (const auto& pair: coeffs) {
            if (pair.first.size() != 2) {
                throw(std::invalid_argument(fmt::format(
                    "Expected quadratic coefficients but got {}\n", 
                    prettyPrint(coeffs)
                )));
            }
        }
        return pure(coeffs); 
    }

    // Trace over the physical part, yielding an exterior element
    //   To compute the trace, for each component convert to majorana basis then compute trace
    RElm tr(const Element& rho) const {
        RElm result = rAlg().zero();
        for (const auto& p: rho.coeffs) {
            auto pivot = p.first.begin() + LAlgebra::num_generators();
            KeyType lkeys(p.first.begin(),pivot), 
                rkeys(pivot, p.first.end());
            auto lelm = car_to_majorana(LElm({{lkeys, FieldType(1.)}}));
            auto accum = RElm({{rkeys, p.second * lelm.tr()}}); 
            result.add_(accum);
        }
        return result; 
    }

    MElm car_to_majorana(const LElm& rho) const {
        auto c = mAlg();
        auto ans = c.zero(); 
        for (auto const& pair:rho.coeffs){
            auto accum = c.one() * pair.second;
            for (uint i=0; i<n; i++) {
                auto a = (c(2*i) + c(2*i+1) * FieldType(0., 1.)) / 2; 
                if (pair.first[2*i] == 1) {
                    accum = accum * a; 
                } 
                if (pair.first[2*i+1] == 1) {
                    accum = accum * a.conj(); 
                }
            }
            ans.add_(accum);
        }
        return ans;
    }

    LElm majorana_to_car(const MElm& rho) const {
        auto a = lAlg();
        auto ans = a.zero(); 
        for (auto const& pair:rho.coeffs){
            auto accum = a.one() * pair.second;
            for (uint i=0; i<n; i++) {
                auto q = (a(2*i) + a(2*i).conj()); 
                auto p = (a(2*i) - a(2*i).conj()) / FieldType(0., 1.);
                if (pair.first[2*i] == 1) {
                    accum = accum * q; 
                } 
                if (pair.first[2*i+1] == 1) {
                    accum = accum * p; 
                }
            }
            ans.add_(accum);
        }
        return ans; 
    }
};

#endif