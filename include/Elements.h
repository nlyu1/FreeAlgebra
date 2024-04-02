#ifndef ELEMENT_H
#define ELEMENT_H
// #include <Eigen/Dense>
#include <fmt/core.h> 
#include "utils.h"
using namespace std; 

// Base class for algebra element operations
template<typename AlgebraRelation>
class AlgebraElement {
public:
    inline static constexpr double SERIES_TOLERANCE = 1e-10;
    inline static constexpr double EQ_TOLERANCE = 1e-7;
    inline static constexpr uint SERIES_MAXITERS = 200;

    // So that we only ever need to instantiate the relation once 
    //    for automorphism relations instantiation is nontrivial 
    static AlgebraRelation& Rel() {
        static AlgebraRelation Rel_{};
        return Rel_;
    }

    CoeffMap coeffs;
    // Constructor taking an rvalue reference (move constructor)
    AlgebraElement(CoeffMap&& map) : 
            coeffs(std::move(map)), n(AlgebraRelation::num_generators()) {
        filter_coeffs_();
        validateKeys();
    }

    AlgebraElement(const CoeffMap& map) : 
            coeffs(map), n(AlgebraRelation::num_generators()) {
        filter_coeffs_();
        validateKeys();
    }

    inline uint num_generators() const { return n; }

    inline static void assert_eqsize(size_t a, size_t b) {
        if (a != b) {
            throw std::invalid_argument(
                fmt::format("Element sizes {} and {} do not match", a, b));
        }
    }

    AlgebraElement clone() const {
        return AlgebraElement(clone_map(coeffs));
    }

    // Construct a zero element (additive identity) of the algebra
    AlgebraElement zero() const {
        return AlgebraElement(CoeffMap({}));
    }

    AlgebraElement one() const {
        return AlgebraElement(CoeffMap({{KeyType(n, 0u), FieldType(1., 0.)}}));
    }
    
    // Scalar operations 
    void add_(const FieldType& scalar) {
        KeyType zero(n, 0u); 
        if (coeffs.contains(zero)) {
            coeffs.at(zero) = coeffs.at(zero) + scalar; 
        } else {
            coeffs[zero] = scalar;
        }
        filter_coeffs_();
    }

    void mul_(const FieldType& scalar) {
        for (const auto& pair:coeffs) {
            coeffs.at(pair.first) = pair.second * scalar;
        }
        filter_coeffs_();
    }

    // Comparison with scalar is considered comparison 
    //      with scalar * multiplicative identity 
    bool operator==(const FieldType& scalar) const {
        return operator==(one() * scalar);
    }

    // Assumes both well-filtered
    bool operator==(const AlgebraElement& other) const {
        return operator-(other).norm() < EQ_TOLERANCE;
    }

    // An element is real if it is its own complex conjugate 
    bool isreal() const {
        return operator==(conj());
    }

    AlgebraElement operator+(const FieldType& scalar) const {
        auto a = clone();
        a.add_(scalar);
        return a; 
    }

    AlgebraElement operator-(const FieldType& scalar) const {
        auto a = clone();
        a.add_(-scalar);
        return a; 
    }

    AlgebraElement operator*(const FieldType& scalar) const {
        auto a = clone();
        a.mul_(scalar);
        return a; 
    }

    AlgebraElement operator/(const FieldType& scalar) const {
        return operator*(FieldType(1, 0) / scalar);
    }

    AlgebraElement operator+(double scalar) const {
        return operator+(FieldType(scalar, 0.0));
    }

    AlgebraElement operator-(double scalar) const {
        return operator-(FieldType(scalar, 0.0));
    }

    AlgebraElement operator*(double scalar) const {
        return operator*(FieldType(scalar, 0.0));
    }

    AlgebraElement operator/(double scalar) const {
        return operator/(FieldType(scalar, 0.0));
    }

    // Pointwise linear addition 
    void add_(const AlgebraElement& other) {
        assert_eqsize(this->n, other.num_generators());
        for (const auto& pair : other.coeffs) {
            if (coeffs.contains(pair.first)) {
                coeffs.at(pair.first) = coeffs.at(pair.first) + pair.second; 
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
        return operator+(other * FieldType(-1, 0));
    }

    AlgebraElement operator*(const AlgebraElement& other) const {
        AlgebraElement result({});
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
        AlgebraElement result({});
        for (const auto& pair: coeffs) {
            KeyType I; 
            auto K = power_to_gen_repr(pair.first); 
            if (K.size() == 0) {
                result.add_(std::conj(pair.second));
                continue; 
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

    // Miscellaneous helpers:
    // Filters null coefficients
    void filter_coeffs_() {
        for (auto it = coeffs.begin(); it != coeffs.end(); ) {
            if (it->second == FieldType(0.0, 0.0)) { // Check if the value is 0
                it = coeffs.erase(it); // Remove the entry and update the iterator
            } else {
                ++it; // Move to the next entry
            }
        }
    }

    std::string to_string() const {
        std::stringstream stream;
        stream << "[\n";
        for (const auto& pair: coeffs) {
            stream << "    ";
            auto I = power_to_gen_repr(pair.first);
            stream << "(";
            for (size_t i = 0; i < I.size(); ++i) {
                if (i > 0) {
                    stream << ", "; 
                }
                stream << Rel().to_string(I[i]);
            }
            stream << "): " << prettyPrint(pair.second) << " \n";
        }
        auto result = stream.str();
        if (coeffs.size() != 0) {
            result.erase(result.length() - 2);
        }
        return result + "\n]";
    }

    // Treating as a vector, collect the absolute norm squared
    double norm() const {
        double result = 0;
        for (const auto& pair: coeffs) {
            result += pair.second.absq();
        }
        return std::pow(result, 0.5);
    }

    FieldType tr() const {
        FieldType result = FieldType(0., 0.);
        for (const auto& pair: coeffs) { 
            result += pair.second * Rel().monomial_tr(pair.first);
        }
        return result;
    }

    AlgebraElement exp() const {
        auto delta = one(), result=one();
        for (uint j=1; j<SERIES_MAXITERS; j++) {
            delta = operator*(delta) / j; 
            result.add_(delta);
            if (j>3 && delta.norm() / result.norm() < SERIES_TOLERANCE) break;
        }
        return result; 
    }

    AlgebraElement log() const {
        auto b = operator-(1.); // b = self - 1
        auto power = one(), result=zero();
        auto this_norm = norm();
        for (uint j=1; j<SERIES_MAXITERS; j++) {
            power = power * b;
            result.add_(power*(std::pow(-1, j+1)/j));
            auto result_norm = result.norm();
            if (j>3 && power.norm()/j / result_norm < SERIES_TOLERANCE) break;
            if (j>3 && result_norm > this_norm * 1e5) {
                throw(std::invalid_argument(
                    std::string("Logarithm does not converge. ") 
                    + std::string("This power series logarithm is not ")
                    + std::string("generally the inverse of exp(). ")
                ));
            }
        }
        return result; 
    }

private:
    uint n; 
    void validateKeys() {
        for (const auto& pair : coeffs) {
            if (pair.first.size() != n) {
                throw std::invalid_argument(
                    fmt::format(
                        "Expected keys to have size {} but got {}", 
                        n, prettyPrint(pair.first)
                    )
                );
            }
        }
    }

    // Figure out the homogeneous exponent relationship from the multiplication relation
    //    (gidx)^pow = FieldType * (gidx)^uint
    tuple<uint, FieldType> homogeneous_exponent(uint gidx, uint pow) const {
        // Handle base cases directly
        if (pow < 2) {
            return {pow, FieldType(1.)}; // Using list initialization for clarity
        }

        // Analyze commutation for homogeneous behavior
        auto t = Rel().commute(gidx, gidx);
        // Square annihilates to 0
        if (t.size() == 0) {
            return {pow, FieldType(0.)}; 
        }
        auto [key, scale] = *(t.begin());
        if (t.size() > 1 || key.size() > 2) {
            throw std::invalid_argument(fmt::format(
                "Expected ({})^2 to be homogeneous in {} but got {}.\n",
                gidx, gidx, prettyPrint(t)));
        }
        FieldType pow_scale = scale.pow(pow - 1);
        uint sq_pow = key.size(), power = 0;

        // Simplify the power based on specific cases
        switch (sq_pow) {
            case 0: power = 0; break; // Generator self-annihilates
            case 1: power = pow % 2; break; // Square of generator annihilates
            case 2: power = pow; break; // No simplification
            default: 
                throw std::invalid_argument(fmt::format(
                    "({})^2 cannot be ({})^{}.\n",
                    gidx, gidx, sq_pow));
        }
        return {power, pow_scale};
    }

    // Given an accumulator and multi-indices in generator multiplication 
    //   representation (and an existing scale), adds to accum the 
    //   corresponding multiplication argument 
    void reorder(const KeyType& I, FieldType scale, 
        AlgebraElement& accum) const {
        if (I.size() == 0) {
            accum.add_(scale);
            return;
        }
        // cout << "Reorder input: " 
        //     << prettyPrint(I) << ": " << scale 
        //     << "    " << accum << endl;
        auto idx = order_violate_idx(I);
        // If in canonical order, arrange homogeneous power then accumulate 
        if (idx == I.size()) {
            FieldType scale_accum = scale; 
            auto J = gen_to_power_repr(I, n);
            // cout << "Original: " << prettyPrint(J) << endl;
            for (size_t i=0; i<J.size(); i++) {
                // Compute the exponent of the J[i]-th power of generator i
                auto t = homogeneous_exponent(i, J[i]);
                scale_accum = scale_accum * std::get<1>(t);
                J[i] = std::get<0>(t); 
                // Shortcut: no change if the scale is already 0 
                if ((scale_accum).abs() == 0.) {
                    return;
                }
            }
            AlgebraElement entry({{J, scale_accum}});
            // cout << "New index: " << prettyPrint(J) << endl;
            // cout << "    Reorder adding:" << entry << endl;
            accum.add_(entry);
            return;
        }
        // AlgebraRelation R;
        // Use commutation relation on the first pair of 
        //  non-canonical product, then recursively call 
        auto L = Rel().commute(I[idx], I[idx + 1]);
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
std::ostream& operator<<(std::ostream& os, const AlgebraElement<T>& element) {
    os << element.to_string(); 
    return os;
}
#endif