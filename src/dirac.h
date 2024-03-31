#ifndef DIRAC_H
#define DIRAC_H
#include "element.h"

class DiracGrassmannElement : public AlgebraElement {
public:
    DiracGrassmannElement(const Eigen::VectorXcd& vec) : AlgebraElement(vec) {}
    DiracGrassmannElement(const std::vector<Complex>& vec) : AlgebraElement(vec) {}

    // Override operator* for pairwise multiplication
    DiracGrassmannElement operator*(const DiracGrassmannElement& other) const {
        // Check if sizes match
        assert_eqsize(this->v.size(), other.v.size());

        // Perform pairwise multiplication
        Eigen::VectorXcd result = this->v.array() * other.v.array();

        return DiracGrassmannElement(result);
    }

    DiracGrassmannElement operator+(const Complex& scalar) const {
        std::cout << type(this->v) << " " << type(std::move(this->v)) << std::endl;
        return DiracGrassmannElement((AlgebraElement(this->v) + scalar).v);
    }
};

#endif