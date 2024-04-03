#ifndef FERMIONS_H
#define FERMIONS_H
#include <cmath>
#include <algorithm>
#include "BaseAlgebraRelations.h"
#include "FinitePowerAlgebras.h"
#include "utils.h"
using namespace std; 


// Majorana algebra on n modes
template<uint n>
class MajoranaAlgebra: public CliffordAlgebra<2*n> {
public:
    using CurAlg = CliffordAlgebra<2*n>;
    using Element = CurAlg::Element; 
    using CARAlg = CARAlgebra<2*n>;
    using CARElement = CARAlg::Element; 
    using SqAlg = ProductAlgebra<CurAlg, CurAlg, decltype(PROD_ANTICOMMUTE)>;
    using SqElement = SqAlg::Element; 

    static CurAlg& curAlg() {
        static CurAlg curAlg_{};
        return curAlg_;
    }

    static SqAlg& sqAlg() {
        static SqAlg sqAlg_{};
        return sqAlg_;
    }

    static CARAlg& carAlg() {
        static CARAlg carAlg_{};
        return carAlg_;
    }

    // Given real coefficients for generators, ensures that 
    //    the result is Hermitian 
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
                coeffs[key] = pair.second;
            } else {
                coeffs[key] = FieldType(0, 1.) * pair.second;
            }
        }
        auto result = Element(coeffs);
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
        return expH / expH.tr(); 
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
        auto U = unitary(coeffs); 
        return U * zero_state() * U.conj();
    }

    // Returns the zero computational (pure) state
    Element zero_state() const {
        auto rho = curAlg().zero();
        auto scale = FieldType(2.).pow(-static_cast<double>(n));
        for (const auto& v: binstr(n)) {
            FieldType term_scale(1.);
            auto term = curAlg().one();
            for (uint i=0; i<n; i++) {
                if (v[i] == 1) {
                    // Each XY=-iZ
                    term_scale = term_scale * FieldType(0., -1);
                    term = term * curAlg()(2*i) * curAlg()(2*i+1);
                }
            }
            rho = rho + term * term_scale * scale; 
        }
        return rho; 
    }

    double degree_weight(const Element& rho, uint k) const {
        double result = 0; 
        for (const auto& pair: rho.coeffs) {
            if (power_to_gen_repr(pair.first).size() == k) {
                result = result + pair.second.absq();
            }
        }
        return result; 
    }

    CoeffMap moments(const Element& rho) {
        CoeffMap coeffs; 
        for (const auto& pair: rho.coeffs) {
            auto deg = power_to_gen_repr(pair.first).size();
            coeffs[pair.first] = pair.second * std::pow(2., n) 
                    * std::pow(-1., deg * (deg - 1) / 2);
        }
        return coeffs; 
    }

    SqElement rotation_unitary(double theta=M_PI/4.) const {
        auto gamma = sqAlg();
        auto H = gamma.zero(); 
        for (uint j=0; j<2*n; j++) {
            H = H + gamma(j) * gamma(2*n + j) * (theta * .5);
        }
        return H.exp();
    }

    Element conv(const Element& rho, const Element& sigma) const {
        auto s = sqAlg();
        auto U = rotation_unitary();
        auto tau = U * s.kron(rho, sigma) * U.conj();
        // cout << "Trace output " << tau << endl;
        return s.trR(tau); 
    }

    CARElement toCAR(const Element& rho) const {
        CARAlg d = carAlg();
        auto result = d.zero();
        for (auto const& p: rho) {
            ;
        }
    }
};

#endif