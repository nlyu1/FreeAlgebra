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
        omega = torch::view_as_complex(omega) * std::pow(2, -.5);

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

    // Cache the noncanonical commutation relation
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

        std::vector<std::vector<FieldType>> nMat; 
        std::vector<std::vector<FieldType>> cMat; 

        /// Initialize nMat: it maps each noncanonical multiplication index 
        //    to its algebraic element (in A's components)
        for (uint i=0; i<n; i++) {
            for (uint j=0; j<=i; j++) {
                // Register the new noncanonical tuple
                uint idx = nspan_inv.size();
                nspan[{i, j}] = idx;
                nspan_inv.push_back({i, j});
                nMat.push_back(std::vector<FieldType>(aspan_inv.size(), FieldType(0., 0.)));

                // Compute the product of b(i) * b(j) in terms of their results in a
                // This is in generator-power representation 
                auto c = (b(i) * b(j)).coeffs;
                for (const auto& pair: c) {
                    auto key = pair.first; 
                    auto value = pair.second; 
                    if (!aspan.contains(key)) {
                        maxdeg = std::max(maxdeg, 
                            static_cast<uint>(power_to_gen_repr(key).size()));
                        // Key is not contained: add the key to all previous entries
                        aspan[key] = aspan.size();
                        aspan_inv.push_back(key);
                        for (auto& pair_: nMat) {
                            pair_.push_back(FieldType(0., 0.));
                        }
                    }
                    // Update the corresponding index 
                    nMat[idx][aspan.at(key)] = value; 
                }
            }
        }
        /// Initialize cmat: it initializes all possible canonical multiplications 
        //    to their canonical A-components. Since the transformation law 
        //    preserves degree, we only need to look for <=maxdeg
        cspan_inv = enumerate_degree_leq(n, maxdeg);
        cout << "nMat initialized. Initializing " << cspan_inv.size() << " candidates." << endl;
        for (uint i=0; i<cspan_inv.size(); i++) {
            // Update the index for the index
            cspan[cspan_inv[i]] = i;
            cMat.push_back(std::vector<FieldType>(aspan_inv.size(), FieldType(0., 0.)));

            auto c = b(cspan_inv[i]).coeffs;
            for (const auto& pair: c) {
                auto key = pair.first; 
                auto value = pair.second; 
                if (!aspan.contains(key)) {
                    // Key is not contained: add key (column) 
                    //    to all previous entries, both for nMat and cMat 
                    aspan[key] = aspan.size();
                    aspan_inv.push_back(key);
                    for (auto& pair_: nMat) {
                        pair_.push_back(FieldType(0., 0.));
                    }
                    for (auto& pair_: cMat) {
                        pair_.push_back(FieldType(0., 0.));
                    }
                }
                // Update the corresponding index 
                cMat[i][aspan.at(key)] = value; 
            }
        }
        cout << "cMat initialized" << endl;
        auto nMat_ = toEigenMatrixXcd(nMat), cMat_overcomplete = toEigenMatrixXcd(cMat);
        cout << "transferred to EigenMatrix" << endl;
        cMat_overcomplete.transposeInPlace();
        nMat_.transposeInPlace();
        // cout << "cMat shape: " << cMat_overcomplete.rows() << ", " << cMat_overcomplete.cols() << endl;
        // Indices corresponding to independent canonical ordering for B. 
        auto indep = independent_cols(cMat_overcomplete); 
        // cout << "Found independent" << endl;
        // cout << cMat_overcomplete << endl;
        std::vector<KeyType> cIdx; // independent canonical indices
        Eigen::MatrixXcd cIndep(cMat_overcomplete.rows(), indep.size());
        for (size_t i = 0; i < indep.size(); i++) {
            cIndep.col(i) = cMat_overcomplete.col(indep[i]);
            cIdx.push_back(cspan_inv[indep[i]]); 
        }
        // cout << "LU decomposition done " << nspan_inv.size() << ", " << cIdx.size() << endl; 
        // Finally! This encodes the noncanonical multiplication solution
        auto multcoeffs = (cIndep.inverse() * nMat_).eval();
        for (uint j=0; j<nspan_inv.size(); j++) { // auto const& nmult: nspan_inv) {
            CoeffMap coeff;
            for (uint i=0; i<cIdx.size(); i++) {
                auto v = FieldType(std::real(multcoeffs(i, j)), std::imag(multcoeffs(i, j))); 
                if (v.absq() < 1e-12) continue; // Ignore numerical instabilities
                coeff[cIdx[i]] = v; 
            }
            noncanonical_commutation[nspan_inv[j][0]][nspan_inv[j][1]] = coeff;
        }
        // cout << "Computation ended." << endl; 

        // cout << "Independent columns " << endl;
        // for (auto i:indep) {
        //     cout << i << " ";
        // }
        // cout << cIndep.rows() << ", " << cIndep.cols() << endl;
        // cout << nMat_.rows() << ", " << nMat_.cols() << endl;
        // cout << multcoeffs << endl;

        // // Debugging: prints out the noncanonical product in A-indices
        // std::stringstream os; 
        // cout << "Independent generators: [";
        // for (auto I: cIdx) {
        //     cout << prettyPrint(power_to_gen_repr(I)) << ", ";
        // }
        // os << "]\naSpan Keys: [";
        // for (uint j0=0; j0<cMat[0].size(); j0++) {
        //     os << prettyPrint(aspan_inv[j0]) << ", ";
        // }
        // os << "]\n";
        // os << "Noncanonical products: [";
        // for (uint j0=0; j0<nspan_inv.size(); j0++) {
        //     os << prettyPrint(nspan_inv[j0]) << ", ";
        // }
        // os << "]\n\n";
        // for (uint i0=0; i0<cMat.size(); i0++) {
        //     auto row = cMat[i0]; 
        //     os << "Canonical B-product: " << prettyPrint(cspan_inv[i0]) << ": ["; 
        //     for (uint j0=0; j0<row.size(); j0++) {
        //         os << prettyPrint(row[j0]) << ", ";
        //     }
        //     os << "]\n";
        // }
        // for (uint i0=0; i0<nMat.size(); i0++) {
        //     auto row = nMat[i0]; 
        //     os << "Noncanonical b-product: " << prettyPrint(nspan_inv[i0]) << ": ["; 
        //     for (uint j0=0; j0<row.size(); j0++) {
        //         os << prettyPrint(row[j0]) << ", ";
        //     }
        //     os << "]\n";
        // }
        // cout << os.str() << endl;

        // for (uint i=0; i<n; i++) {
        //     for (uint j=0; j<=i; j++) {
        //         cout << i << " " << j << ": " 
        //         << prettyPrint(power_to_gen_repr(noncanonical_commutation[i][j])) << endl;
        //     }
        // }
    }

    CoeffMap commute_noncanonical(uint i, uint j) const override {
        return noncanonical_commutation[i][j]; 
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
            ans.add_(a(j) * FieldType(J()[idx][j]));
        }
        return ans; 
    }

    // Returns the A-components of a multi-index in b-coordinates
    ElmA b(KeyType I) {
        auto ans = ElmA(CoeffMap()).one();
        for (uint i=0; i<n; i++) {
            if (ans.norm() < 1e-10) {
                return ans;
            }
            ans = ans * b(i).pow(I[i]); 
        }
        return ans;
    }
};

#endif