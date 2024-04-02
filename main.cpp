// #include "FinitePowerAlgebras.h"
// #include "Automorphism.h"
// #include "Fermions.h"
// #include <torch/torch.h>
// #include<iostream>

// using namespace std;

// template<typename GammaType>
// auto z(const GammaType& gamma) {
//     return gamma(0).conj() * gamma(1) + gamma(3);
// }

// typedef CliffordAlgebra<3> LAlg;
// typedef CARAlgebra<4> RAlg; 

// int main(){
//     // LAlg alpha;
//     // RAlg beta;
//     // ProductAlgebra<LAlg, RAlg, decltype(PROD_ANTICOMMUTE)> gamma; 
//     // cout << alpha(0) << endl;
//     // cout << beta(0) << endl;
//     // cout << gamma.extR(alpha(0)) << gamma.extL(beta(0)) << endl;
//     // cout << gamma.anticommutator(gamma.extR(alpha(0)), gamma.extL(beta(0))) << endl; 
//     // cout << (gamma.projL(gamma.extR(alpha(0)) + gamma.extL(beta(0)))) << endl;
//     // cout << gamma.extL(beta(0)) << endl;
//     // cout << (gamma.projL(gamma.extL(beta(0)))) << endl;
//     // cout << (gamma.projR(gamma.extR(alpha(0))) == beta.zero()) << endl;

//     const uint n=3;
//     MajoranaAlgebra<n> gamma;
//     // ProductAlgebra<MajoranaAlgebra<n>, MajoranaAlgebra<n>, decltype(PROD_ANTICOMMUTE)> gamma_sq; 

//     RealCoeffMap coeffs = {
//         {{}, 3.},
//         {{0, 1}, 1.}, 
//         {{1, 2}, 2.},
//         {{2, 3}, 3.}
//     };

//     // Gaussian hamiltonian 
//     auto rho = gamma.density(coeffs);
//     auto U = gamma.unitary(coeffs);
//     // cout << "New   " << endl;
//     // cout << gamma.moments(rho) << endl;
//     // cout << rho << endl;
//     // cout << gamma_sq.kron(rho, rho) << endl;
//     // cout << gamma_sq.kron(rho, rho).tr() << endl;

//     // const uint n=3; 
//     // CoeffMap coeffs = {
//     //     {{}, Complex(3.)},
//     //     {{0, 1}, Complex(1.)}, 
//     //     {{1, 2}, Complex(2.)},
    
//     //     {{2, 3}, Complex(3.)}
//     // };
//     // MajoranaAlgebra<n> gamma; 
//     // cout << gamma(0).conj() << endl;

//     DiracMajoranaBijection<6> R;
//     CARAlgebra<6> a;
//     ImageRelation<CanonicalAnticommRelation<6>, DiracMajoranaBijection<6>> r;
// }



#include "utils.h"
#include "BaseAlgebraRelations.h"
#include "FinitePowerAlgebras.h"
#include "Automorphism.h"
#include<iostream>
using namespace std; 

// // Test algebraic relations 
// int main() {
//     ImageRelation<CanonicalAnticommRelation<4>, DiracMajoranaBijection<4>> a;



// }

int main() {
    // Create a 2x2 complex tensor
    auto real = torch::tensor({{1., 2.}, 
                                {3., 4.}});
    auto imag = torch::tensor({{5., 6.}, 
                                {7., 8.}});
    auto omega = torch::cat({real.unsqueeze(-1), imag.unsqueeze(-1)}, -1);
    omega = torch::view_as_complex(omega);
    omega.set_requires_grad(true);

    omega.sum().backward();
    cout << omega.grad() << endl;

    // // Perform operations
    // auto element1 = tensor.index({0, 0}); // Extracting first element
    // auto element2 = tensor.index({1, 1}); // Extracting second element
    // auto product = element1 * element2;   // Multiplying the elements

    // // Perform differentiation with respect to the product
    // product.backward();

    // // Inspect the gradient of the original tensor
    // std::cout << "Gradient of the tensor: " << tensor.grad() << std::endl;

    return 0;
}