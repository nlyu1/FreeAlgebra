#ifndef FERMIONS_H
#define FERMIONS_H
#include <cmath>
#include "BaseAlgebraRelations.h"
#include "FinitePowerAlgebras.h"
using namespace std; 

// Majorana algebra on n modes
template<uint n>
class MajoranaAlgebra: public CliffordAlgebra<2*n> {
public:
    using Element = CliffordAlgebra<2*n>::Element; 

    // All except the identity is traceless. 
    // Identity has trace 2**(n/2)
    Complex tr(const Element& X) const {
        KeyType zero(2*n, 0u);
        return X.coeffs.contains(zero) ? X.coeffs.at(zero) * std::pow(2., n) : 0.;
    }

    // Given real coefficients for generators, ensures that 
    //    the result is Hermitian 
    Element hamiltonian(const RealCoeffMap& coeffs_) const {
        CoeffMap coeffs;
        for (const auto& pair: coeffs_) {
            auto degree = pair.first.size();
            if (degree * (degree - 1) / 2 % 2 == 0){ 
                coeffs[gen_to_power_repr(pair.first, 2*n)] = pair.second;
            } else {
                coeffs[gen_to_power_repr(pair.first, 2*n)] = pair.second * Complex(0, 1.);
            }
        }
        auto result = Element(coeffs);
        cout << result << endl;
        assert (result.isreal()); 
        return result;
    }

    // Assembles a unitary from the coefficient of its exponents 
    Element unitary(const RealCoeffMap& coeffs) const {
        return (hamiltonian(coeffs) * Complex(0., 1.)).exp();
    }

    // Assembles a thermal state from its Hamiltonian 
    Element density(const RealCoeffMap& coeffs) const {
        auto expH = hamiltonian(coeffs).exp();
        return expH / tr(expH); 
    }

    double degree_weight(const Element& rho, uint k) const {
        double result = 0; 
        for (const auto& pair: rho.coeffs) {
            if (power_to_gen_repr(pair.first).size() == k) {
                result += std::pow(std::abs(pair.second)*pow(2., n), 2);
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

    CliffordAlgebra<2*n>::Element rotation_unitary(double theta=M_PI/4.) const {
        CliffordAlgebra<2*n> gamma; 
        auto H = gamma.zero();
        for (uint j=0; j<2*n; j++) {
            H = H + gamma(j) * gamma(n + j) * (theta * .5);
        }
        return H.exp();
    }

    Element conv(const Element& rho, const Element& sigma) const {
        ;
    }
};

#endif