#ifndef AUTOMORPHISM_H
#define AUTOMORPHISM_H
// #include "Elements.h"
#include "BaseAlgebraRelations.h"
#include <torch/torch.h>
using namespace std;

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

        uint j = 2, i = 1;
        cout << prettyPrint(SrcRel().commute_noncanonical(j, i)) << endl;
    }
    
};


#endif