#ifndef UTILS_H
#define UTILS_H

#include <complex>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <tuple>

typedef std::complex<double> Complex;
typedef uint PowType;
typedef std::vector<PowType> KeyType;
typedef Complex ValueType;
struct VectorHash {
    std::size_t operator()(const KeyType& v) const {
        std::size_t hash = 0;
        for (PowType i : v) {
            hash ^= std::hash<PowType>{}(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

typedef std::unordered_map<KeyType, ValueType, VectorHash> CoeffMap; 


// Pretty printing of a double value
std::string prettyPrint(double value) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(6) << value; // Set a max precision
    std::string str = stream.str();
    str.erase(str.find_last_not_of('0') + 1, std::string::npos);
    // If the decimal point is now the last character, remove it as well
    if (str.back() == '.') {
        str.erase(str.length() - 1);
    }
    return str;
}

// String formatting of complex number
std::string prettyPrint(Complex z) {
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
    return result.empty() ? "0" : "("+result+")";
}

bool isPowerOfTwo(size_t n) {
        return n > 0 && (n & (n - 1)) == 0;
}

// Binary and integer conversion 
uint32_t binaryToInt(const std::vector<int>& binary) {
    uint32_t result = 0;
    int power = 0;
    for (int i = binary.size() - 1; i >= 0; --i) {
        result += binary[i] * static_cast<uint32_t>(pow(2, power));
        ++power;
    }
    return result;
}

std::vector<int> intToBinary(uint32_t k, int n) {
    std::vector<int> binary(n, 0);
    for (int i = n - 1; i >= 0 && k > 0; --i) {
        binary[i] = k % 2;
        k /= 2;
    }
    return binary;
}

std::string prettyPrint(const KeyType& x) {
    std::stringstream stream;
    stream << "[";
    for (const auto& v: x) {
        stream << v << ","; 
    }
    auto result = stream.str();
    if (x.size() > 0) {
        result.erase(result.length() - 1);
    }
    return result + "]";
}


std::string prettyPrint(const CoeffMap& map){
    std::stringstream stream;
    stream << "{";
    for (const auto& pair: map) {
        stream 
        << "(" << prettyPrint(pair.first)
        << ", " << prettyPrint(pair.second) << "), ";
    }
    auto result = stream.str();
    if (map.size() != 0) {
        result.erase(result.length() - 2);
    }
    return result + "}";
}

CoeffMap value_map(const CoeffMap& map, 
    std::function<ValueType(const ValueType&)> func) {
    CoeffMap newMap; 
    for (const auto& pair:map) {
        newMap[pair.first] = func(pair.second);
    }
    return newMap; 
}

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
#endif
