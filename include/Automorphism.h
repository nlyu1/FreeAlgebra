#ifndef AUTOMORPHISM_H
#define AUTOMORPHISM_H
// #include "Elements.h"
#include "BaseAlgebraRelations.h"
#include <torch/torch.h>
using namespace std;

typedef std::unordered_map<KeyType, uint, VectorHash> FlexibleIdx;

template<uint n>
class DiracMajoranaBijection {
public:
    static_assert(n % 2 == 0, 
        "Dirac to Majorana bijection requires an even number of generators");
    torch::Tensor J, JInv; 

    DiracMajoranaBijection<n>() {
        auto real = torch::tensor({{1., 1.}, 
                                   {0., 0.}});
        auto imag = torch::tensor({{0., 0.}, 
                                   {-1., 1.}});
        // Use torch::cat to concatenate real and imag parts along a new dimension
        auto omega = torch::cat({real.unsqueeze(-1), imag.unsqueeze(-1)}, -1);
        omega = torch::view_as_complex(omega);

        J = torch::block_diag(std::vector<torch::Tensor>(n / 2, omega));
        JInv = J.inverse();
    }

    static constexpr uint num_generators() {
        return n;
    }
};


template<typename SrcRelation, typename Bijection>
struct ImageRelation : BaseRelation<SrcRelation::num_generators()> {
    static_assert(SrcRelation::num_generators() == Bijection::num_generators(), 
        "Source algebra and bijection-rule must specify the same number of generators");
    static constexpr auto n = SrcRelation::num_generators();

    static SrcRelation& SrcRel() {
        static SrcRelation Rel{};
        return Rel;
    }
    static torch::Tensor& J() {
        static Bijection B{};
        return B.J; 
    }
    using ElmA = AlgebraElement<SrcRelation>;

    // Caches the noncanonical commutation relation
    CoeffMap noncanonical_commutation[n][n];

    ImageRelation() {
        // Compute the canonical anticommutation relation for the image algebra. 
        //    Linear bijection preserves degree, so b_i*b_j = C_{ijkl} a_ka_k
        //    Set up 2n linear equations
        auto real = torch::randn({2*n, 2*n});
        auto imag = torch::randn({2*n, 2*n});
        auto E = torch::complex(torch::randn({2*n, 2*n}), \
            torch::randn({2*n, 2*n}));  // Equation coefficients
        // (matrix of noncanonical ordering) = (conversion matrix) * (matrix of canonical ordering)

        uint maxdeg = 2; // Maximum degree that a noncanonical product can blow up to
        // Mapping from span of the noncanonical products of b to those of a
        // Records the span of the noncanonical products in terms of a's generators
        FlexibleIdx nspan; 
        std::vector<KeyType> nspan_inv; 

        // Index for the degrees of A involved in everything
        FlexibleIdx aspan; 
        std::vector<KeyType> aspan_inv; 
        
        // Mapping from span of the noncanonical products of b to those of a
        FlexibleIdx cspan; 
        std::vector<KeyType> cspan_inv; 

        std::vector<std::vector<ValueType>> nMat; 
        std::vector<std::vector<ValueType>> cMat; 
        for (uint i=0; i<n; i++) {
            for (uint j=0; j<=i; j++) {
                // Register the new noncanonical tuple
                uint idx = nspan_inv.size();
                nspan[{i, j}] = idx;
                nspan_inv.push_back({i, j});
                nMat.push_back(std::vector<ValueType>(aspan_inv.size(), Complex(0., 0.)));

                // Compute the product of b(i) * b(j) in terms of their results in a


                for (const auto& pair: SrcRel().commute(i, j)) {
                    auto key = pair.first; 
                    auto value = pair.second; 
                    // If key contained: just modify the entry 
                    if (!nspan.contains(key)) {
                        // If key is not contained: add the key 
                        aspan[key] = aspan.size();
                        aspan_inv.push_back(key);
                        for (auto& pair_: nMat) {
                            pair_.push_back(Complex(0., 0.));
                        }
                    }
                    // Update the corresponding index 
                    nMat[idx][aspan.at(key)] = value; 
                }
                // cout << fmt::format("({}, {}) @ {}: \n", i, j, idx) 
                // << prettyPrint(nMat) << endl;
            }
        }

        std::stringstream os; 
        os << "Keys: [";
        for (uint j0=0; j0<nMat[0].size(); j0++) {
            os << prettyPrint(aspan_inv[j0]) << ", ";
        }
        os << "]\n";
        for (uint i0=0; i0<nMat.size(); i0++) {
            auto row = nMat[i0]; 
            os << "Noncanonical product: " << prettyPrint(nspan_inv[i0]) << ": ["; 
            for (uint j0=0; j0<row.size(); j0++) {
                os << prettyPrint(row[j0]) << ", ";
            }
            os << "]\n";
        }
        cout << os.str() << endl;

        cout << b(0) << endl;
    }

private:
    ElmA a(uint i) const {
        assert(i >= 0 && i< SrcRelation::num_generators()); 
        KeyType coeffs({i});
        return ElmA({{gen_to_power_repr(coeffs, SrcRelation::num_generators()), 1}});
    }

    ElmA b(uint idx) {
        auto ans = ElmA(CoeffMap());
        for (uint j=0; j<n; j++) {
            auto c = J()[idx][j];
            // double real_part = c.real().item<double>();
            // double imag_part = c.imag().item<double>();
            // ans.add_(a(j) * Complex(c.real().item<double>, c.imag().item<double>));
        }
        return ans; 
    }
};

#endif