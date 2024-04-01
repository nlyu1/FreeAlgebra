#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include <pybind11/complex.h>
#include <pybind11/operators.h>
#include <unordered_map>
#include "Elements.h"

namespace py = pybind11;

PYBIND11_MODULE(freealgebra, m) {
    // py::class_<AlgebraElement>(m, "AlgebraElement");
    //     // Miscellaneous operations
    //     .def("n", &AlgebraElement::n)
    //     .def("v", &AlgebraElement::coeffs)
    //     // Scalar operations
    //     .def("__add__", static_cast<AlgebraElement (AlgebraElement::*)(const Complex&) const>(&AlgebraElement::operator+))
    //     .def("__sub__", static_cast<AlgebraElement (AlgebraElement::*)(const Complex&) const>(&AlgebraElement::operator-))
    //     .def("__mul__", static_cast<AlgebraElement (AlgebraElement::*)(const Complex&) const>(&AlgebraElement::operator*))
    //     .def("__truediv__", static_cast<AlgebraElement (AlgebraElement::*)(const Complex&) const>(&AlgebraElement::operator/))
    //     // Vector operations
    //     .def("__add__", static_cast<AlgebraElement (AlgebraElement::*)(const AlgebraElement&) const>(&AlgebraElement::operator+))
    //     .def("__sub__", static_cast<AlgebraElement (AlgebraElement::*)(const AlgebraElement&) const>(&AlgebraElement::operator-));

    // py::class_<DiracGrassmannElement, AlgebraElement>(m, "DiracGrassmannElement")
    //     .def(py::init<const std::vector<Complex>&>())
    //     .def("__mul__", &DiracGrassmannElement::operator*);

        // .def("n", &DiracGrassmannElement::n)
        // .def("v", &DiracGrassmannElement::coeffs)

        // Scalar operations
        // .def("__add__", static_cast<DiracGrassmannElement (DiracGrassmannElement::*)(const Complex&) const>(&AlgebraElement::operator+))
        // .def("__sub__", static_cast<DiracGrassmannElement (DiracGrassmannElement::*)(const Complex&) const>(&DiracGrassmannElement::operator-))
        // .def("__mul__", static_cast<DiracGrassmannElement (DiracGrassmannElement::*)(const Complex&) const>(&DiracGrassmannElement::operator*))
        // .def("__truediv__", static_cast<DiracGrassmannElement (DiracGrassmannElement::*)(const Complex&) const>(&DiracGrassmannElement::operator/))
        
        // // Vector operations
        // .def("__add__", static_cast<DiracGrassmannElement (DiracGrassmannElement::*)(const DiracGrassmannElement&) const>(&DiracGrassmannElement::operator+))
        // .def("__sub__", static_cast<DiracGrassmannElement (DiracGrassmannElement::*)(const DiracGrassmannElement&) const>(&DiracGrassmannElement::operator-));
}