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


#endif
