#ifndef UTILS_H
#define UTILS_H

#include <complex>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <vector>
#include <tuple>
#include <Eigen/Dense>
#include <fmt/core.h> 
#include <torch/torch.h>
#include <unordered_map>
#include "complex.h"
using namespace std; 


typedef uint PowType;
typedef std::vector<PowType> KeyType;
struct VectorHash {
    std::size_t operator()(const KeyType& v) const {
        std::size_t hash = 0;
        for (PowType i : v) {
            hash ^= std::hash<PowType>{}(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

typedef std::unordered_map<KeyType, FieldType, VectorHash> CoeffMap; 
typedef std::unordered_map<KeyType, double, VectorHash> RealCoeffMap;


// String formatting of complex number
std::string prettyPrint(FieldType z) {
    std::string result;
    double re = z.real();
    double im = z.imag();
    if (re != 0.0) {
        result += prettyPrint(re);
    }
    if (im != 0.0) {
        if (!result.empty() && im > 0) {
            result += "+";
        } else if (im < 0) {
            result += "-";
            im = -im; // Make positive for printing
        }
        result += prettyPrint(im) + "i";
    }
    if (result.empty()) {
        return result+"0";
    } else if (re == 0. || im == 0.) {
        return result; 
    } else {
        return "("+result+")";
    }
}

std::string prettyPrint(const KeyType& x) {
    std::stringstream stream;
    stream << "[";
    for (const auto& v: x) {
        stream << prettyPrint(v) << ","; 
    }
    auto result = stream.str();
    if (x.size() > 0) {
        result.erase(result.length() - 1);
    }
    return result + "]";
}


template<typename T>
std::string prettyPrint(const std::unordered_map<KeyType, T, VectorHash>& coeffs) {
        std::stringstream stream;
        stream << "[\n";
        for (const auto& pair: coeffs) {
            stream << "    ";
            auto I = pair.first;
            stream << "(";
            for (size_t i = 0; i < I.size(); ++i) {
                if (i > 0) {
                    stream << ", "; 
                }
                stream << I[i];
            }
            stream << "): " << prettyPrint(pair.second) << " \n";
        }
        auto result = stream.str();
        if (coeffs.size() != 0) {
            result.erase(result.length() - 2);
        }
        return result + "\n]\n";
    }

// template<typename T>
// std::string prettyPrint(const std::vector<std::vector<T>>& matrix) {
//     std::stringstream os; 
//     for (const auto& row : matrix) {
//         for (const auto& element : row) {
//             os << prettyPrint(element) << ' ';
//         }
//         os << fmt::format("[{}]\n", row.size()); // End each row with a new line
//     }
//     return os.str(); // Return the ostream object to chain the operator calls
// }


CoeffMap clone_map(const CoeffMap& map) {
    CoeffMap newMap; 
    for (const auto& pair: map) {
        newMap[pair.first] = pair.second;
    }
    return newMap;
}

// Returns the first index which is greater than its subsequent element 
uint order_violate_idx(const KeyType& K) {
    for (uint i=0; i<K.size()-1; i++) {
        if (K[i] > K[i+1]) {
            return i;
        }
    }
    return K.size();
}
// Checks that a generator-multiplication representation is well-ordered
void assert_well_ordered(const KeyType& K) {
    assert (order_violate_idx(K)==K.size());
}

// Converts generator multiplication representation 
// to the power representation
KeyType gen_to_power_repr(const KeyType& K, uint n) {
    KeyType ans(n, 0);
    if (K.size() == 0) {
        return ans; 
    }
    assert_well_ordered(K);
    for (const auto& k: K) {
        assert(k >= 0 && k < n);
        ans[k] += 1;
    }
    return ans; 
}

KeyType power_to_gen_repr(const KeyType& K) {
    KeyType ans;
    uint idx = 0;
    for (const auto& k: K) {
        if (k != 0) {
            for (uint i=0; i<k; i++) {
                ans.push_back(idx);
            }
        }
        idx += 1;
    }
    return ans; 
}

template<typename T>
std::vector<T> mergeVectors(const std::vector<T>& I, const std::vector<T>& J) {
    if (I.empty()) {
        return J; 
    } else if (J.empty()) {
        return I; 
    } else {
        std::vector<T> merged = I;
        merged.insert(merged.end(), J.begin(), J.end());
        return merged;
    }
}

std::ostream& operator<<(std::ostream& os, const CoeffMap& coeffs) {
    os << prettyPrint(coeffs);
    return os;
}

template <typename T>
bool allzero(std::vector<T> v) {
    for (auto const& x:v) {
        if (x != 0) {
            return false; 
        }
    }
    return true; 
}

// Returns a vector of vectors which returns all n-tuples of sum deg
std::vector<KeyType> enumerate_degree_eq(uint n, uint deg) {
    if (n == 1) {
        return {{deg}}; 
    } else {
        std::vector<KeyType> result; 
        for (uint d=0; d<=deg; d++) {
            auto rec = enumerate_degree_eq(n-1, d);
            for (uint k=0; k<rec.size(); k++) {
                rec[k].push_back(deg - d);
            }
            result.insert(result.end(), rec.begin(), rec.end());
        }
        return result; 
    }
}

std::vector<KeyType> enumerate_degree_leq(uint n, uint deg) {
    std::vector<KeyType> result;
    for (uint m=0; m<=deg; m++){
        auto rec = enumerate_degree_eq(n, m);
        result.insert(result.end(), rec.begin(), rec.end());
    }
    return result;
}


Eigen::MatrixXcd toEigenMatrixXcd(const vector<vector<FieldType>>& matrix) {
    const size_t rows = matrix.size();
    const size_t cols = matrix.empty() ? 0 : matrix[0].size();
    // Initialize an Eigen matrix of complex doubles
    Eigen::MatrixXcd eigenMatrix(rows, cols);
    // Fill the matrix
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            eigenMatrix(i, j) = matrix[i][j].item(); // Directly assign the complex number
        }
    }
    return eigenMatrix;
}

KeyType independent_cols(const Eigen::MatrixXcd& matrix) {
    Eigen::FullPivLU<Eigen::MatrixXcd> lu_decomp(matrix);
    Eigen::VectorXi indices = lu_decomp.permutationQ().indices();

    KeyType independentCols;
    for(int i = 0; i < lu_decomp.rank(); ++i) {
        independentCols.push_back(static_cast<uint>(indices(i)));
    }
    return independentCols;
}

#endif
