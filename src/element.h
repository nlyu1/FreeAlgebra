#ifndef ELEMENT_H
#define ELEMENT_H
// #include <vector>
// #include <complex>
#include <Eigen/Dense>
#include <fmt/core.h> 
#include "utils.h"
using namespace std; 

template<typename CommutationRelation>
class AlgebraElement {
public:
    CoeffMap coeffs;
    uint n; // Number of free generators

    // Constructor taking an rvalue reference (move constructor)
    AlgebraElement(CoeffMap&& map, uint n)
        : coeffs(std::move(map)), n(n)  {
        validateKeys();
    }

    // Constructor taking a const lvalue reference
    AlgebraElement(const CoeffMap& map, uint n)
        : coeffs(map), n(n) {
        validateKeys();
    }

    std::string tostring() {
        return prettyPrint(this->coeffs);
    }

    inline uint num_generators() { return n; }

    inline static void assert_eqsize(size_t a, size_t b) {
        if (a != b) {
            throw std::invalid_argument(fmt::format("Element sizes {} and {} do not match", a, b));
        }
    }

    AlgebraElement clone() const {
        return AlgebraElement(clone_map(coeffs), n);
    }
    
    // Scalar operations 
    void add_(const Complex& scalar) {
        KeyType zero(n, 0u); 
        if (coeffs.contains(zero)) {
            coeffs.at(zero) += scalar; 
        } else {
            coeffs[zero] = scalar;
        }
    }

    void mul_(const Complex& scalar) {
        for (const auto& pair:coeffs) {
            coeffs.at(pair.first) = pair.second * scalar;
        }
    }

    AlgebraElement operator+(const Complex& scalar) const {
        auto a = clone();
        a.add_(scalar);
        return a; 
    }

    AlgebraElement operator-(const Complex& scalar) const {
        auto a = clone();
        a.add_(-scalar);
        return a; 
    }

    AlgebraElement operator*(const Complex& scalar) const {
        auto a = clone();
        a.mul_(scalar);
        return a; 
    }

    AlgebraElement operator/(const Complex& scalar) const {
        return operator*(Complex(1, 0) / scalar);
    }

    // Pointwise linear addition 
    void add_(const AlgebraElement& other) {
        assert_eqsize(this->n, other.n);
        for (const auto& pair : other.coeffs) {
            if (coeffs.contains(pair.first)) {
                coeffs.at(pair.first) += pair.second; 
            } else {
                coeffs[pair.first] = pair.second; 
            }
        }
    }

    AlgebraElement operator+(const AlgebraElement& other) const {
        auto a = clone();
        a.add_(other);
        return a;
    }

    AlgebraElement operator-(const AlgebraElement& other) const {
        return operator+(other * Complex(-1, 0));
    }

    AlgebraElement operator*(const AlgebraElement& other) const {
        AlgebraElement result({}, n);
        for (const auto& pair1: coeffs) {
            for (const auto& pair2: other.coeffs) {
                reorder(
                    mergeVectors(
                        power_to_gen_repr(pair1.first), 
                        power_to_gen_repr(pair2.first)
                    ), pair1.second * pair2.second, result);
            }
        }
        return result;
    }

private:
    void validateKeys() {
        for (const auto& pair : coeffs) {
            if (pair.first.size() != n) {
                throw std::invalid_argument(
                    fmt::format(
                        "Expected keys to have size {} but got {}", n, prettyPrint(pair.first)
                    )
                );
            }
        }
    }
    // Given an accumulator and multi-indices in generator multiplication 
    //   representation (and an existing scale), adds to accum the 
    //   corresponding multiplication argument 
    void reorder(const KeyType& I, ValueType scale, 
        AlgebraElement& accum) const {
        if (I.size() == 0) {
            accum.add_(scale);
            return;
        }
        // cout << "Reorder input: " 
        //     << prettyPrint(I) << ": " << scale 
        //     << "    " << prettyPrint(accum) << endl;

        auto idx = order_violate_idx(I);
        // If in canonical order, then add to accum
        if (idx == I.size()) {
            AlgebraElement entry({{gen_to_power_repr(I, n), scale}}, n);
            // cout << "    Reorder adding:" << prettyPrint(entry) << endl;
            accum.add_(entry);
            return;
        }
        CommutationRelation cr;
        // Use commutation relation on the first pair of 
        //  non-canonical product, then recursively call 
        auto L = cr.commute(I[idx], I[idx + 1]);
        for (const auto& pair: L) {
            auto newI = KeyType(I.begin(), I.begin()+idx);
            newI.insert(newI.end(), pair.first.begin(), pair.first.end());
            newI.insert(newI.end(), I.begin()+idx+2, I.end());
            reorder(newI, scale*pair.second, accum);
        }
        return;
    }
};

template<typename T>
std::string prettyPrint(const AlgebraElement<T>& a) {
    return prettyPrint(a.coeffs);
}


// Grassmann commutation relation 
struct GCR {
    CoeffMap commute(uint i, uint j) const {
        CoeffMap result;
        if (i < j) {
            throw(std::invalid_argument("Expects non-canonical ordering"));
        }
        if (i > j) {
            result[{j, i}] = Complex(-1., 0.);
        }
        return result; 
    }
};

// Dirac commutation relation
struct DCR {
    CoeffMap commute(uint i, uint j) const {
        CoeffMap result;
        if (i < j) {
            throw(std::invalid_argument("Expects non-canonical ordering"));
        }
        if (i > j) {
            if (i % 2 == 1 && j == i - 1) {
                result[{}] = Complex(1., 0.);
            }   
            result[{j, i}] = Complex(-1., 0.);
        }
        return result; 
    }
};

typedef AlgebraElement<DCR> DiracElm;
typedef AlgebraElement<GCR> ExtElm;
#endif