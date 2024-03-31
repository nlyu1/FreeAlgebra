#ifndef ELEMENT_H
#define ELEMENT_H
#include <Eigen/Dense>
#include <fmt/core.h> 
#include "utils.h"
using namespace std; 

template<typename AlgebraRelation>
class AlgebraElement {
public:
    CoeffMap coeffs;
    uint n; // Number of free generators

    // Constructor taking an rvalue reference (move constructor)
    AlgebraElement(CoeffMap&& map, uint n)
        : coeffs(std::move(map)), n(n)  {
        filter_coeffs_();
        validateKeys();
    }

    // Constructor taking a const lvalue reference
    AlgebraElement(const CoeffMap& map, uint n)
        : coeffs(map), n(n) {
        filter_coeffs_();
        validateKeys();
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

    // Construct a zero element (additive identity) of the algebra
    AlgebraElement zero() const {
        return AlgebraElement(CoeffMap({}), n);
    }

    AlgebraElement one() const {
        return AlgebraElement(CoeffMap({{KeyType(n, 0u), Complex(1., 0.)}}), n);
    }
    
    // Scalar operations 
    void add_(const Complex& scalar) {
        KeyType zero(n, 0u); 
        if (coeffs.contains(zero)) {
            coeffs.at(zero) += scalar; 
        } else {
            coeffs[zero] = scalar;
        }
        filter_coeffs_();
    }

    void mul_(const Complex& scalar) {
        for (const auto& pair:coeffs) {
            coeffs.at(pair.first) = pair.second * scalar;
        }
        filter_coeffs_();
    }

    // Comparison with scalar is considered comparison with scalar * multiplicative identity 
    bool operator==(const Complex& scalar) const {
        return operator==(one() * scalar);
    }

    // Assumes both well-filtered
    bool operator==(const AlgebraElement& other) const {
        return coeffs == other.coeffs;
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
        filter_coeffs_();
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

    AlgebraElement pow(uint k) {
        auto result = one(); // Multiplicative identity
        // cout << "Self:" << prettyPrint(self) << endl;
        for (uint j=0; j<k; j++){
            result = operator*(result);
        }
        return result;
    }

    // Conjugation relation
    AlgebraElement conj() const {
        AlgebraRelation R;
        AlgebraElement result({}, n);
        for (const auto& pair: coeffs) {
            KeyType I; 
            auto K = power_to_gen_repr(pair.first); 
            if (K.size() == 0) {
                result.add_(std::conj(pair.second));
                return result; 
            }
            for (int j=K.size()-1; j>=0; j--){
                I.push_back(R.conj(K[j]));
            }
            // cout << "Reordering: " << prettyPrint(I) << endl;
            reorder(I, std::conj(pair.second), result);
        }
        result.filter_coeffs_();
        return result; 
    }

    // Filters null coefficients
    void filter_coeffs_() {
        for (auto it = coeffs.begin(); it != coeffs.end(); ) {
            if (it->second == Complex(0.0, 0.0)) { // Check if the value is 0
                it = coeffs.erase(it); // Remove the entry and update the iterator
            } else {
                ++it; // Move to the next entry
            }
        }
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
            // Enforce canonical anticommutation relation now: duplicates go to 0
            //  In the future: enforce non-fixed points
            AlgebraElement entry({{gen_to_power_repr(I, n), scale}}, n);
            for (size_t i=0; i<I.size() - 1; i++){ // Zero duplicates
                if (I[i] == I[i+1]) {
                    entry = entry * Complex(0., 0.);
                }
            }
            // cout << "    Reorder adding:" << prettyPrint(entry) << endl;
            accum.add_(entry);
            return;
        }
        AlgebraRelation R;
        // Use commutation relation on the first pair of 
        //  non-canonical product, then recursively call 
        auto L = R.commute(I[idx], I[idx + 1]);
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
#endif