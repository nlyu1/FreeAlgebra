#include "utils.h"
#include "BaseAlgebraRelations.h"
#include "FinitePowerAlgebras.h"
#include "Automorphism.h"
#include "ProductPowerAlgebra.h"
#include "Fermions.h"
#include<iostream>
#include <cassert>
#include<cmath>
using namespace std; 

const uint n = 5;

double h(uint a) {
    return (a>n/2)?-1.:1.;
}

double h(uint a, uint b) {
    if (a == b) {
        return 0.;
    } else if (a > b) {
        return -h(b, a);
    } else {
        return 2*a+b;
    }
}

double g(uint a, uint b) {
    if (a == b) {
        return 0.;
    } else if (a > b) {
        return -g(b, a);
    } else {
        return (a/2==0 && b/3==0)? b*2+a : 0.;
    }
}

// Constant displacement 
double g(uint k) {
    return k;
}

void print_h() {
    cout << "Matrix: "; 
    for (auto i=0; i<n; i++) {
        cout << endl;
        for (auto j=0; j<n; j++) {
            cout << h(i, j) << " ";
        }
    }
    cout << endl;

    cout << "Displacement ";
    for (auto j=0; j<n; j++) {
        cout << h(j) << " ";
    }
    cout << endl;
}

void print_g() {
    cout << "Matrix: "; 
    for (auto i=0; i<n; i++) {
        cout << endl;
        for (auto j=0; j<n; j++) {
            cout << g(i, j) << " ";
        }
    }
    cout << endl;

    cout << "Displacement ";
    for (auto j=0; j<n; j++) {
        cout << g(j) << " ";
    }
    cout << endl;
}

double delta(uint a, uint b) {
    return (a==b)? 1. : 0.; 
}

// int main() {
//     CliffordAlgebra<n> gamma;
//     print_h();

//     cout << "G" << endl;
//     print_g();

//     auto hquad = gamma.zero();
//     for (uint a=0; a<n; a++) {
//         for (uint b=0; b<n; b++) {
//             hquad = hquad + gamma(a) * gamma(b) * (h(a, b) / 2);
//         }
//     }

//     auto gquad = gamma.zero();
//     for (uint a=0; a<n; a++) {
//         for (uint b=0; b<n; b++) {
//             gquad = gquad + gamma(a) * gamma(b) * (g(a, b) / 2);
//         }
//     }

//     auto glin = gamma.zero();
//     for (uint a=0; a<n; a++) {
//         glin = glin + gamma(a) * g(a);
//     }

//     auto hlin = gamma.zero();
//     for (uint a=0; a<n; a++) {
//         hlin = hlin + gamma(a) * h(a);
//     }
//     hlin = hlin * FieldType(0., -1.); // Times -i

//     // Now, onto verification of the commutators 
//     auto qlcom_true = gamma.commutator(hquad, glin); 
//     auto qlcom = gamma.zero();
//     for (uint k=0; k<n; k++) {
//         for (uint a=0; a<n; a++) {
//             qlcom = qlcom - gamma(a) * (h(k, a) * g(k) * 2); 
//         }
//     }
//     cout << "Quadratic-linear commutator consistency: " << (qlcom == qlcom_true) << endl;

//     auto lqcom_true = gamma.commutator(hlin, gquad); 
//     auto lqcom = gamma.zero();
//     for (uint k=0; k<n; k++) {
//         for (uint a=0; a<n; a++) {
//             lqcom = lqcom + gamma(a) * (FieldType(0., -2.) * g(k, a) * h(k)); 
//         }
//     }
//     cout << "Quadratic-linear commutator consistency: " << (lqcom == lqcom_true) << endl;

//     auto llcom_true = gamma.commutator(hlin, glin); 
//     auto llcom = gamma.zero(); 
//     for (uint j=0; j<n; j++) {
//         for (uint k=0; k<n; k++) {
//             if (j == k) continue; 
//             llcom = llcom + gamma(j) * gamma(k) * (FieldType(0., 2.) * h(k) * g(j));
//         }
//     }
//     cout << "Linear-linear commutator consistency: " << (llcom == llcom_true) << endl;


//     auto qqcom_true = gamma.commutator(hquad, gquad); 
//     auto qqcom = gamma.zero(); 

//     for (uint a=0; a<n; a++) {
//         for (uint b=0; b<n; b++) {
//             for (uint u=0; u<n; u++) {
//                 if (u == b) continue;
//                 auto bracket = 
//                     // gamma(b) * gamma(v) * delta(a, u) * (1. - delta(b, v)) * h(b, a) * g(u, v) + 
//                     gamma(b) * gamma(u) * h(a, b) * g(u, a) * 2;
//                 qqcom = qqcom + bracket;
//             }
//         }
//     }
//     cout << "Quadratic-quadratic commutator consistency: " << (qqcom == qqcom_true) << endl;
// }



int main() {
    CliffordAlgebra<n> gamma;
    print_h();

    cout << "G" << endl;
    print_g();

    auto hquad = gamma.zero();
    for (uint a=0; a<n; a++) {
        for (uint b=0; b<n; b++) {
            hquad = hquad + gamma(a) * gamma(b) * (h(a, b) / 2);
        }
    }
    auto gquad = gamma.zero();
    for (uint a=0; a<n; a++) {
        for (uint b=0; b<n; b++) {
            gquad = gquad + gamma(a) * gamma(b) * (g(a, b) / 2);
        }
    }
    auto glin = gamma.zero();
    for (uint a=0; a<n; a++) {
        glin = glin + gamma(a) * g(a);
    }
    auto hlin = gamma.zero();
    for (uint a=0; a<n; a++) {
        hlin = hlin + gamma(a) * h(a);
    }
    hlin = hlin * FieldType(0., -1.); // Times -i

    // Now, onto verification of the commutators 
    auto qlcom_true = gamma.commutator(hquad, glin); 
    auto qlcom = gamma.zero();
    for (uint k=0; k<n; k++) {
        for (uint a=0; a<n; a++) {
            qlcom = qlcom - gamma(a) * (h(k, a) * g(k) * 2); 
        }
    }
    cout << "Quadratic-linear commutator consistency: " << (qlcom == qlcom_true) << endl;

    auto lqcom_true = gamma.commutator(hlin, gquad); 
    auto lqcom = gamma.zero();
    for (uint k=0; k<n; k++) {
        for (uint a=0; a<n; a++) {
            lqcom = lqcom + gamma(a) * (FieldType(0., -2.) * g(k, a) * h(k)); 
        }
    }
    // cout << lqcom_true << endl << lqcom << endl;
    cout << "Quadratic-linear commutator consistency: " << (lqcom == lqcom_true) << endl;

    auto llcom_true = gamma.commutator(hlin, glin); 
    auto llcom = gamma.zero(); 
    for (uint j=0; j<n; j++) {
        for (uint k=0; k<n; k++) {
            if (j == k) continue; 
            llcom = llcom + gamma(j) * gamma(k) * (FieldType(0., 2.) * h(k) * g(j));
        }
    }
    // cout << llcom_true << endl;
    // cout << llcom << endl;
    cout << "Linear-linear commutator consistency: " << (llcom == llcom_true) << endl;


    auto qqcom_true = gamma.commutator(hquad, gquad); 
    auto qqcom = gamma.zero(); 
    // for (uint a=0; a<n; a++) {
    //     for (uint b=0; b<n; b++) {
    //         for (uint u=0; u<n; u++) {
    //             for (uint v=0; v<n; v++) {
    //                 // auto bracket = 
    //                 //     gamma(b) * gamma(v) * -2. * delta(a, u) * (1. - delta(b, v)) + 
    //                 //     gamma(a) * gamma(v) * 2. * delta(b, u) * (1. - delta(a, v)) + 
    //                 //     gamma(b) * gamma(u) * 2 * delta(a, v) * (1. - delta(b, u)) + 
    //                 //     gamma(a) * gamma(u) * -2. * delta(b, v) * (1. - delta(a, u)); 
    //                 // qqcom = qqcom + bracket * (h(a, b) * g(u, v) * .25);
    //                 // Second-level equivalence 
    //                 auto bracket = 
    //                     // gamma(b) * gamma(v) * delta(a, u) * (1. - delta(b, v)) * h(b, a) * g(u, v) + 
    //                     gamma(b) * gamma(u) * delta(a, v) * (1. - delta(b, u)) * h(a, b) * g(u, v) * 2;
    //                 qqcom = qqcom + bracket;
    //                 // auto true_bracket = gamma.commutator(gamma(a) * gamma(b), gamma(u) * gamma(v)); 
    //                 // No need to care about these cases, since the coefficients are zero 
    //                 // if ((a == b) || (u == v)) {
    //                 //     continue; 
    //                 // }
    //                 // if (bracket != true_bracket) {
    //                 //     cout << fmt::format("({}, {}, {}, {})", a, b, u, v) << bracket << true_bracket << endl;
    //                 // }
    //             }
    //         }
    //     }
    // }

    for (uint a=0; a<n; a++) {
        for (uint b=0; b<n; b++) {
            for (uint u=0; u<n; u++) {
                if (u == b) continue;
                auto bracket = 
                    // gamma(b) * gamma(v) * delta(a, u) * (1. - delta(b, v)) * h(b, a) * g(u, v) + 
                    gamma(b) * gamma(u) * h(a, b) * g(u, a) * 2;
                qqcom = qqcom + bracket;
            }
        }
    }
    // cout << qqcom_true << endl;
    // cout << qqcom << endl;
    cout << "Quadratic-quadratic commutator consistency: " << (qqcom == qqcom_true) << endl;
}


// int main() {
//     const uint n = 2;
//     DiracAlgebra<n> a; 
//     auto beta_coeffs = ComplexDispType({
//         FieldType(1., 2.), 
//         FieldType(0., 1.)
//     }), gamma_coeffs = ComplexDispType({
//         FieldType(1., 0.), FieldType(.5, .3)
//     });
//     auto beta_vec = a.disp_vec(beta_coeffs), gamma_vec = a.disp_vec(gamma_coeffs); 
//     auto beta_cohr = a.coherent(beta_coeffs); 
//     auto gamma_cohr = a.coherent(gamma_coeffs); 
//     auto overlap_sq = a.tr(beta_cohr * gamma_cohr);
//     auto Dbeta = a.disp(beta_coeffs), Dgamma = a.disp(gamma_coeffs);  
//     cout << "Computed overlap_sq: " << overlap_sq << endl;
//     cout << "Trace equality: " 
//         << (a.tr(beta_cohr * gamma_cohr) == a.tr(gamma_cohr * beta_cohr)) << endl;
//     cout << "Coherence check: " 
//         << (beta_cohr == Dbeta * a.vac() * Dbeta.conj()) << ", " 
//         << (gamma_cohr == Dgamma * a.vac() * Dgamma.conj()) << endl;
//     cout << "Overlap expansion check: " << 
//         (overlap_sq == a.tr(Dbeta * a.vac() * Dbeta.conj() * Dgamma * a.vac() * Dgamma.conj())) << endl;
//     cout << "Overlap expansion cyclic check: " << 
//         (overlap_sq == a.tr(Dbeta * a.vac() * Dbeta.conj() * Dgamma * a.vac() * Dgamma.conj())) << endl;
//     // Predicted overlap: 
//     auto pred_overlap_sq = a.zero();
//     for (uint i=0; i<n; i++) {
//         pred_overlap_sq.add_((beta_vec[i] - gamma_vec[i]) 
//             * (beta_vec[i].conj() - gamma_vec[i].conj())); 
//     }
//     pred_overlap_sq = pred_overlap_sq.exp(); 
//     cout << "Predicted overlap_sq: " << pred_overlap_sq << endl;
// }

// int main() {
//     const uint n = 4;
//     DiracAlgebra<n> a; 
//     auto beta_coeffs = RealCoeffMap({
//         {{0, 1, 4, 5}, 1}, 
//         {{0, 2, 4, 6}, 1},
//     }); 
//     auto rho = a.pure(beta_coeffs);
//     cout << rho << rho.isreal() << endl;
//     cout << "Majorana representation " 
//         << a.car_to_majorana(a.projL(rho)) << endl;
//     auto k = a.zero();
//     auto l = a.lAlg(); 
//     auto r = a.rAlg();
//     for (uint i=0; i<n; i++) {
//         auto ai = a.extR(l(2*i)); 
//         auto xii = a.extL(r(2*i)); 
//         // Glauber 113
//         k.add_(xii * ai.conj() + xii.conj() * ai);
//     }
//     cout << "Computed k: " << k << endl;
//     // k = k.exp(); 
//     auto d = a.disp(ComplexDispType(n, FieldType(1.))); 
//     cout << "Commutator: " << a.commutator(rho, k) << endl;
//     auto xi = a.tr(rho * k);
//     // cout << "CAR exterior: " << xi << endl;
//     // cout << "Maj exterior: " << a.car_to_majorana_ext(xi) << endl;
//     // cout << "Consistency: " 
//     //     << (xi == a.majorana_to_car_ext(a.car_to_majorana_ext(xi))) << endl;
//     cout << "Xi: " << xi << endl;
//     auto zeta = a.log(xi); 
//     cout << "Cumulants: " << zeta << zeta.isreal() << endl;
//     auto w = a.fourier_transform(xi); 
//     cout << "Wigner function: " << w << w.isreal() << endl;
//     auto logw = a.log(w); 
//     cout << "Log-wigner function: " << logw << logw.isreal() << endl;
//     auto entropy = (w * logw); 
//     cout << "Entropy: " << entropy << endl;
// }