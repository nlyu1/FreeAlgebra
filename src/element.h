#ifndef ELEMENT_H
#define ELEMENT_H
// #include <vector>
// #include <complex>
#include <Eigen/Dense>
#include <fmt/core.h> 
#include "utils.h"

class AlgebraElement {
public:
    Eigen::VectorXcd v; // VectorXcd is a vector of complex doubles in Eigen

    // Internal copy constructor from Eigen-array
    AlgebraElement(const Eigen::VectorXcd& vec) : v(vec) {}
    AlgebraElement(Eigen::VectorXcd&& vec) : v(std::move(vec)) {}

    // Python-facing constructor from std::vector
    AlgebraElement(const std::vector<Complex>& vec) : v(vec.size()) {
        if (!isPowerOfTwo(vec.size())) {
            throw std::invalid_argument(
                fmt::format(
                    "Expected AlgebraElement coefficient to be power of 2 but received length {}"
                    , vec.size())
            );
        }
        for (size_t i = 0; i < vec.size(); ++i) {
            v(i) = vec[i];
        }
    }

    uint32_t n() { return static_cast<uint32_t>(std::log2(v.size()));}

    std::vector<Complex> coeffs() {
        std::vector<Complex> result; 
        for (auto j = 0; j < v.size(); j++) {
            result.push_back(v(j)); 
        }
        return result; 
    }

    // Vector operations 
    static void assert_eqsize(size_t a, size_t b) {
        if (a != b) {
            throw std::invalid_argument(fmt::format("Element sizes {} and {} do not match", a, b));
        }
    }
    
    // Scalar operations 
    AlgebraElement operator+(const Complex& scalar) const {
        return AlgebraElement(v.array() + scalar);
    }
    AlgebraElement operator-(const Complex& scalar) const {
        return AlgebraElement(v.array() - scalar);
    }
    AlgebraElement operator*(const Complex& scalar) const {
        return AlgebraElement(v.array() * scalar);
    }
    AlgebraElement operator/(const Complex& scalar) const {
        return AlgebraElement(v.array() / scalar);
    }
    // Operator overload for pointwise linear operations between two AlgebraElements
    AlgebraElement operator+(const AlgebraElement& other) const {
        assert_eqsize(this->v.size(), other.v.size());
        return AlgebraElement(v + other.v);
    }
    AlgebraElement operator-(const AlgebraElement& other) const {
        assert_eqsize(this->v.size(), other.v.size());
        return AlgebraElement(v - other.v);
    }
    AlgebraElement operator*(const AlgebraElement& other) const {
        throw std::invalid_argument("Algebraic multiplication undefined");
    }
};
#endif